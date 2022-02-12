#!env python3
# Time-Stamp: <2022-02-14 21:24:25>

# Copyright 2022 Sho Iwamoto / Misho
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#     http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging
import pathlib
import math
import sys
import coloredlogs
import yaslha

# For illustrative purposes.
import tokenize
import re
import sys

from typing import Union, Optional, List, Any, Dict, Tuple, Sequence
import collections.abc

import convert_slha1_to_slha2

logger = logging.getLogger(__name__)

# for type-hinting
ParamList = Dict[str, Optional[float]]

DEFAULT = {
    "YU": [6.8e-6, 0.0341, 0.931],
    "YD": [1.47e-5, 0.000293, 0.01553],
    "YE": [2.7930e-6, 0.00058838, 0.0099944],
}

# Dictionary keys (block names) should be in capital letters
# MATRIX_BLOCKS = {
#     "UPMNS": ("RMNS", 3, 3),
#     "IMUPMNS": ("IMNS", 3, 3),
#     "VCKM": ("RCKM", 3, 3),
#     "IMVCKM": ("ICKM", 3, 3),
#     "NMIX": ("RNN", 4, 4),
#     "IMNMIX": ("INN", 4, 4),
#     "UMIX": ("RUU", 2, 2),
#     "IMUMIX": ("IUU", 2, 2),
#     "VMIX": ("RVV", 2, 2),
#     "IMVMIX": ("IVV", 2, 2),
#     "SNUMIX": ("RRn", 3, 3),
#     "IMSNUMIX": ("IRn", 3, 3),
#     "SELMIX": ("RRl", 6, 6),
#     "IMSELMIX": ("IRl", 6, 6),
#     "USQMIX": ("RRu", 6, 6),
#     "IMUSQMIX": ("IRu", 6, 6),
#     "DSQMIX": ("RRd", 6, 6),
#     "IMDSQMIX": ("IRd", 6, 6),
#     "YU": ("Ryu", 3, 3),
#     "IMYU": ("Iyu", 3, 3),
#     "YD": ("Ryd", 3, 3),
#     "IMYD": ("Iyd", 3, 3),
#     "YE": ("Rye", 3, 3),
#     "IMYE": ("Iye", 3, 3),
#     "MSL2": ("RmL2", 3, 3),
#     "IMMSL2": ("ImL2", 3, 3),
#     "MSE2": ("RmE2", 3, 3),
#     "IMMSE2": ("ImE2", 3, 3),
#     "MSQ2": ("RmQ2", 3, 3),
#     "IMMSQ2": ("ImQ2", 3, 3),
#     "MSU2": ("RmU2", 3, 3),
#     "IMMSU2": ("ImU2", 3, 3),
#     "MSD2": ("RmD2", 3, 3),
#     "IMMSD2": ("ImD2", 3, 3),
#     "TE": ("Rte", 3, 3),
#     "IMTE": ("Ite", 3, 3),
#     "TU": ("Rtu", 3, 3),
#     "IMTU": ("Itu", 3, 3),
#     "TD": ("Rtd", 3, 3),
#     "IMTD": ("Itd", 3, 3),
# }
# PLAIN_BLOCKS = {
#     "FRALPHA": {0: "alp"},  # apply special treatment
#     "SMINPUTS": {1: "aEWM1", 3: "aS"},
#     "HMIX": {1: "RMUH", 2: "tb", 4: "MA2"},
#     "IMHMIX": {1: "IMUH"},
#     "MSOFT": {1: "RMx1", 2: "RMx2", 3: "RMx3", 21: "mHu2", 22: "mHd2"},
#     "IMMSOFT": {1: "IMx1", 2: "IMx2", 3: "IMx3"},
# }  # type: Dict[str, Dict[int, str]]


class Model:
    ignore_blocks = ["SPINFO", "DCINFO", "MODSEL", "MINPAR", "EXTPAR"]

    def __init__(self):
        # type: ()->None
        self.blocks = {}  # type: Dict[str, Dict[Tuple[int, ...], str]]
        self.masses = {}  # type: Dict[int, str]
        self.widths = {}  # type: Dict[int, str]
        self.default_values = {}  # type: Dict[str, float]

    def add_parameter(self, p):  # type: (Parameter)->None
        self.default_values[p.name] = float(p.value)

        if p.lhablock is None:
            logger.error(f"UFO parameter {p.name} does not have LHABLOCK property.")
            exit(1)

        assert p.lhablock.upper() == p.lhablock
        if not isinstance(p.lhacode, collections.abc.Sequence) or not all(
            isinstance(i, int) for i in p.lhacode
        ):
            logger.error(f"UFO parameter {p.name} has LHACODE not a list of INT.")
            exit(1)

        if p.lhablock == "MASS":
            assert len(p.lhacode) == 1
            self.masses[p.lhacode[0]] = p.name
        elif p.lhablock == "DECAY":
            assert len(p.lhacode) == 1
            self.widths[p.lhacode[0]] = p.name
        else:
            key = tuple(i for i in p.lhacode)
            if p.lhablock not in self.blocks:
                self.blocks[p.lhablock] = {}
            assert key not in self.blocks[p.lhablock]
            self.blocks[p.lhablock][key] = p.name


class Parameter:
    target_model = Model()

    def __init__(self, name, nature, type, value, texname, lhablock=None, lhacode=None):
        # type: (str, str, str, Union[float, str], str, Any, Any)->None
        self.name = name
        self.nature = nature
        self.type = type
        self.value = value
        self.texname = texname
        self.lhablock = lhablock.upper() if lhablock else None
        self.lhacode = lhacode
        if self.nature.lower() == "external":
            self.target_model.add_parameter(self)


def load_model(model_dir):  # type: (pathlib.Path)->Model
    parameters_file = model_dir / "parameters.py"
    if not parameters_file.is_file():
        logger.error(f"File not found: {parameters_file}")
        exit(1)
    parameters_definition = parameters_file.read_text()
    parameters_definition = re.sub(r"^from .*$", "", parameters_definition, flags=re.M)

    exec(parameters_definition)
    return Parameter.target_model


def parse_slha(slha, model):  # type: (yaslha.slha.SLHA, Model) -> List[str]
    new_parameters = dict((k, v) for k, v in model.default_values.items())
    is_set = dict((k, False) for k in model.default_values)

    def set(key, value, ignore_missing=False):
        if key in new_parameters:
            new_parameters[key] = value
            is_set[key] = True
        elif ignore_missing:
            return
        else:
            logger.error(
                f"SLHA has a parameter {key}={value} but {key} is not in the UFO."
            )
            logger.info(f"key: {key}, value: {value}")
            exit(1)

    for block in slha.blocks.values():
        block_name = block.name.upper()
        if block_name in model.ignore_blocks:
            pass
        elif block_name == "ALPHA":
            set(model.blocks["FRALPHA"][(1,)], block[None])
        elif block_name == "MASS":
            for key, value in block.items():
                set(model.masses[key], value)
        elif block_name in model.blocks:
            definition = model.blocks[block_name]
            for k, v in block.items():
                if isinstance(k, int):
                    key = (k,)
                else:
                    key = tuple(i for i in k)
                if param_name := definition.get(key):
                    set(param_name, v)
                elif len(key) >= 2 and v == 0.0:
                    # zero elements in matrices are safely ignored.
                    logger.debug(f"Ignored zero: {block_name} {k}")
                else:
                    if block_name == "SMINPUTS":
                        logger.info(f"Ignored: {block_name} {k} {block.comment[k]}")
                    else:
                        logger.warning(f"Ignored: {block_name} {k} {block.comment[k]}")
        else:
            logger.warning(f"Block ignored: {block_name}")
    for pid, decay in slha.decays.items():
        set(model.widths[pid], decay.width)

    # check missing parameters
    missing_messages = {}  # type: Dict[str, str]
    missing_messages_critical = []  # type: List[str]
    for key, value in new_parameters.items():
        if is_set[key]:
            continue
        if m := re.match(r"R(MNS|CKM)(\d)x(\d)", key):
            if m[1] not in missing_messages:
                missing_messages[m[1]] = f"Model-default for {m[1]} is used: "
            missing_messages[m[1]] += f"{m[2]}{m[3]}={value} "
        elif m := re.match(r"R(yu|yd|ye)(\d)x(\d)", key):
            y_name = m[1].upper()
            if y_name not in missing_messages:
                missing_messages[y_name] = f"Converter-default for {y_name} is used: "
            if m[2] == m[3]:
                sm_yukawa = DEFAULT[y_name][int(m[2]) - 1]
                beta = math.atan(new_parameters["tb"])
                beta_factor = math.sin(beta) if y_name == "YU" else math.cos(beta)
                new_parameters[key] = sm_yukawa / beta_factor
            else:
                new_parameters[key] = 0.0
            missing_messages[y_name] += f"{m[2]}{m[3]}={new_parameters[key]} "
        elif m := re.match(r"M(Z|W|e|m|ta|U|C|T|D|S|B)", key):
            if "mass" not in missing_messages:
                missing_messages["mass"] = f"Model-default SM masses are used: "
            missing_messages["mass"] += f"{m[1]}={value} "
        elif m := re.match(r"W(Z|W|H|T)", key):
            if "width" not in missing_messages:
                missing_messages["width"] = f"Model-default SM widths are used: "
            missing_messages["width"] += f"{m[1]}={value} "
        elif key.startswith("W"):
            if value == 0.0:
                missing_messages[key] = f"Default ZERO-WIDTH is used: {key}={value}"
            else:
                if "w" not in missing_messages:
                    missing_messages["w"] = f"Model-default non-zero MSSM width for "
                missing_messages["w"] += key[1:] + " "
        else:
            missing_messages_critical.append(
                f"Model-default parameter may have potential issue: {key}={value}"
            )
    result = [f"{key} = {value}" for key, value in new_parameters.items()]
    for message in missing_messages.values():
        logger.info(message)
        result.append("# Info: " + message)
    for message in missing_messages_critical:
        logger.warning(message)
        result.append("# Warning: " + message)
    return result


if __name__ == "__main__":
    coloredlogs.install(logger=logging.getLogger(), fmt="%(levelname)8s %(message)s")

    slha1_mode = False
    if len(sys.argv) == 4 and (
        sys.argv[1] == "--slha1" or sys.argv[1] == "--slha1-fullmix"
    ):
        slha1_mode = True
        model_dir = pathlib.Path(sys.argv[2])
        slha = yaslha.parse_file(sys.argv[3])
    elif len(sys.argv) == 3:
        model_dir = pathlib.Path(sys.argv[1])
        slha = yaslha.parse_file(sys.argv[2])
    else:
        logger.critical(f"Usage: {sys.argv[0]} (--slha1) model_directory SLHA_file")
        exit(1)
    model = load_model(model_dir)
    if slha1_mode:
        slha = convert_slha1_to_slha2.convert_slha1_to_slha2(
            slha, sys.argv[1] == "--slha1-fullmix"
        )
    result_lines = parse_slha(slha, model)
    for i in result_lines:
        print(i)
