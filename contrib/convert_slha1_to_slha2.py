#!env python3
# Time-Stamp: <2022-02-16 01:08:17>

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

import copy
import logging
import math
import sys
from typing import List, Tuple

import yaslha

logger = logging.getLogger(__name__)

# Default value of Yukawa couplings used for AU-TU conversion.
# These values are in the SM and thus must be divided by sin(beta) or cos(beta).
DEFAULT = {
    "YU": [6.8e-6, 0.0341, 0.931],
    "YD": [1.47e-5, 0.000293, 0.01553],
    "YE": [2.7930e-6, 0.00058838, 0.0099944],
}

KEPT_BLOCKS = [
    "SPINFO",
    "DCINFO",
    "MODSEL",
    "SMINPUTS",
    "MINPAR",
    "EXTPAR",
    "MASS",  # MASS and NMIX need to be fixed for negative neutralino mass
    "NMIX",
    "UMIX",
    "VMIX",
    "HMIX",
    "ALPHA",
    "GAUGE",
    "YU",
    "YD",
    "YE",
]


def set_t(slha1, slha2, t_name, a_name, y_name, beta_factor):
    # type: (yaslha.slha.SLHA, yaslha.slha.SLHA, str, str, str, float)->None
    a_block = slha1[a_name]
    t_block = copy.deepcopy(a_block)
    t_block.head.name = t_name
    slha2.blocks[t_name] = t_block

    for (i, j), value in slha1[a_name].items():
        if not (i == j and i in [1, 2, 3]):
            logger.warning(f"Unaccepted entry {a_name}{i}{j} found; ignored.")
            continue
        if value == 0.0:
            t_block[i, j] = 0.0

        else:
            y = slha1.get(y_name, i, j)
            if y is None:
                y = DEFAULT[y_name][i - 1] / beta_factor
                logger.warning(f"Entry {y_name}{i}{i} missing; use default value {y}.")
            t_block[i, j] = value * y  # see Eq.6 of hep-ph/0311123
            t_block.comment[i, j] = f"updated; A={value}"


def handle_msoft(slha1, slha2, msoft):
    # type: (yaslha.slha.SLHA, yaslha.slha.SLHA, yaslha.block.Block)->None

    def create_block(name):  # type: (str)->None
        if name not in slha2.blocks:
            block = yaslha.block.Block(name)
            block.head.comment = "converted from MSOFT"
            slha2.add_block(block)
            if msoft.q is not None:
                slha2[name].q = msoft.q

    for k, v in msoft.items():
        if k in [1, 2, 3, 21, 22]:
            create_block("MSOFT")
            slha2["MSOFT"][k] = v
            slha2["MSOFT"].comment[k] = msoft.comment[k]
        elif 31 <= k <= 36 or 41 <= k <= 49:
            if k in [31, 32, 33]:
                (target, key) = ("MSL2", k - 30)
            elif k in [34, 35, 36]:
                (target, key) = ("MSE2", k - 33)
            elif k in [41, 42, 43]:
                (target, key) = ("MSQ2", k - 40)
            elif k in [44, 45, 46]:
                (target, key) = ("MSU2", k - 43)
            elif k in [47, 48, 49]:
                (target, key) = ("MSD2", k - 46)
            else:
                # comes here only if k is not an integer.
                raise
            create_block(target)
            slha2[target][key, key] = v * v * 1.0
            slha2[target].comment[key, key] = f"MSOFT {k} = {v}"
        else:
            logger.warning(f"Unaccepted entry MSOFT {k} found; ignored.")


def handle_2x2_mixing(slha1_block, slha2_block, generations):
    # type: (yaslha.blocks.Block, yaslha.blocks.Block, Tuple[int, int])->None

    def is_close(x, y):  # type: (float, float)->bool
        return abs(1.0 - x / y) < 1e-4

    # check validity (although not quite necessary)
    cos1 = cos2 = sin1 = sin2 = 0.0
    for k, v in slha1_block.items():
        if k == (1, 1):
            cos1 = v
        elif k == (1, 2):
            sin1 = v
        elif k == (2, 1):
            sin2 = v
        elif k == (2, 2):
            cos2 = v
        else:
            logger.warning(f"Unaccepted entry {slha1_block.name} {k} found; ignored.")
    if not (
        (cos2 == 0 or is_close(abs(cos1 / cos2), 1.0))
        and (sin2 == 0 or is_close(abs(sin1 / sin2), 1.0))
        and is_close(cos1 ** 2 + sin1 ** 2, 1.0)
        and not (sin1 * sin2 * cos1 * cos2 > 0)
    ):
        logger.warning(f"SLHA1 block {slha1_block.name} is accepted but seems strange.")

    # sfermion sorting (by mass) is done in the later steps.
    for i, i_gen in enumerate(generations):
        for j, j_gen in enumerate(generations):
            slha2_block[i_gen, j_gen] = slha1_block.get(i + 1, j + 1, default=0.0) * 1.0
            slha2_block.comment[i_gen, j_gen] = "converted from " + slha1_block.name


def flip_neutralino_negative_mass(slha2):
    # type: (yaslha.slha.SLHA)->None
    neutralinos = [(1, 1000022), (2, 1000023), (3, 1000025), (4, 1000035)]
    for i, pid in neutralinos:
        if slha2["MASS", pid] < 0:
            slha2["MASS", pid] = abs(slha2["MASS", pid])
            for block_name in ["NMIX", "IMNMIX"]:
                if block_name not in slha2.blocks:
                    new_block = yaslha.block.Block(block_name)
                    for i2 in (1, 2, 3, 4):
                        for j2 in (1, 2, 3, 4):
                            new_block[i2, j2] = 0.0
                    slha2.add_block(new_block)

            for j in (1, 2, 3, 4):
                new_r = (-1) * slha2.get("IMNMIX", i, j, default=0.0)
                new_i = slha2.get("NMIX", i, j, default=0.0)
                slha2["NMIX"][i, j] = new_r
                slha2["IMNMIX"][i, j] = new_i
                slha2["NMIX"].comment[i, j] = slha2["IMNMIX"].comment[i, j] = "flipped"


def reorder_sfermion_mass_matrices(slha2, full_mixing):
    # type: (yaslha.slha.SLHA, bool)->None

    def sub(mixing_name, pids):
        # type: (str, List[int])->None
        length = len(pids)
        p_list = []
        for i, pid in enumerate(pids):
            p_list.append(
                (
                    slha2["MASS", pid],
                    [slha2[mixing_name][i + 1, j + 1] for j in range(length)],
                    [slha2[mixing_name].comment[i + 1, j + 1] for j in range(length)],
                )
            )
        if full_mixing:
            p_list.sort(key=lambda x: x[0])
        elif len(pids) == 6:
            # generation-specific mixing
            for i in range(3):
                if p_list[i][0] > p_list[i + 3][0]:
                    (p_list[i], p_list[i + 3]) = (p_list[i + 3], p_list[i])
        else:
            pass  # nothing for sneutirnos

        for i, pid in enumerate(pids):
            slha2["MASS"][pid] = p_list[i][0]
            for j in range(length):
                slha2[mixing_name][i + 1, j + 1] = p_list[i][1][j]
                slha2[mixing_name].comment[i + 1, j + 1] = p_list[i][2][j]

    sub("USQMIX", [1000002, 1000004, 1000006, 2000002, 2000004, 2000006])
    sub("DSQMIX", [1000001, 1000003, 1000005, 2000001, 2000003, 2000005])
    sub("SELMIX", [1000011, 1000013, 1000015, 2000011, 2000013, 2000015])
    sub("SNUMIX", [1000012, 1000014, 1000016])


def convert_slha1_to_slha2(slha1, full_mixing):
    # type: (yaslha.slha.SLHA, bool)->yaslha.slha.SLHA
    """
    Return SLHA2-formatted object equivalent to SLHA1 input.

    If full_mixing is true, sfermions are sorted by its mass order.
    Otherwise, first, second, and third generation fermions are set in (1,4), (2,5), and
    (3,6)-th entry, respectively.

    Following unofficial SLHA1 blocks are accepted:
     - SUPMIX, SDOWNMIX, SCHARMMIX, SSTRANGEMIX, SEMIX, SMUMIX
    """
    slha2 = yaslha.slha.SLHA()
    beta = math.atan(slha1["HMIX"][2])

    for (gen, name) in [(6, "DSQMIX"), (6, "USQMIX"), (6, "SELMIX"), (3, "SNUMIX")]:
        block = yaslha.block.Block(name)
        block.head.comment = "converted from SLHA1 mixing blocks"
        for i in range(1, gen + 1):
            for j in range(1, gen + 1):
                block[i, j] = 1.0 if i == j else 0.0
                block.comment[i, j] = "(preset)"
        slha2.add_block(block)

    for block_name, block in slha1.blocks.items():
        block_name = block_name.upper()
        if block_name in KEPT_BLOCKS:
            slha2[block_name] = block
        elif block_name == "AU":
            set_t(slha1, slha2, "TU", "AU", "YU", math.sin(beta))
        elif block_name == "AD":
            set_t(slha1, slha2, "TD", "AD", "YD", math.cos(beta))
        elif block_name == "AE":
            set_t(slha1, slha2, "TE", "AE", "YE", math.cos(beta))
        elif block_name == "MSOFT":
            handle_msoft(slha1, slha2, block)
        elif block_name == "SUPMIX":
            handle_2x2_mixing(block, slha2["USQMIX"], (1, 4))
        elif block_name == "SCHARMMIX":
            handle_2x2_mixing(block, slha2["USQMIX"], (2, 5))
        elif block_name == "STOPMIX":
            handle_2x2_mixing(block, slha2["USQMIX"], (3, 6))
        elif block_name == "SDOWNMIX":
            handle_2x2_mixing(block, slha2["DSQMIX"], (1, 4))
        elif block_name == "SSTRANGEMIX":
            handle_2x2_mixing(block, slha2["DSQMIX"], (2, 5))
        elif block_name == "SBOTMIX":
            handle_2x2_mixing(block, slha2["DSQMIX"], (3, 6))
        elif block_name == "SEMIX":
            handle_2x2_mixing(block, slha2["SELMIX"], (1, 4))
        elif block_name == "SMUMIX":
            handle_2x2_mixing(block, slha2["SELMIX"], (2, 5))
        elif block_name == "STAUMIX":
            handle_2x2_mixing(block, slha2["SELMIX"], (3, 6))
        else:
            logger.warning(f"Unknown SLHA1 block {block_name} found; copied to SLHA2")
            slha2[block_name] = block

    flip_neutralino_negative_mass(slha2)
    reorder_sfermion_mass_matrices(slha2, full_mixing)

    return slha2


if __name__ == "__main__":
    full_mixing = False
    file = None
    if len(sys.argv) == 3 and sys.argv[1] == "--full":
        full_mixing = True
        file = sys.argv[2]
    elif len(sys.argv) == 2:
        file = sys.argv[1]
    if file is None:
        logger.critical(f"Usage: {sys.argv[0]} [--full] input_path")
        exit(1)
    slha = yaslha.parse_file(file)
    slha2 = convert_slha1_to_slha2(slha, full_mixing)
    dumper = yaslha.dumper.SLHADumper()
    dumper.set_config("comments_preserve", yaslha.dumper.CommentsPreserve.TAIL)
    print(dumper.dump(slha2))
