# @Date:   2019-10-29T17:38:25+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: TagePDB.py
# @Last modified time: 2019-10-29T17:42:48+08:00
import pandas as pd
import numpy as np
import sys
from MMCIFplus import MMCIF2DictPlus
sys.path.append("../py_src/")

_NECESSARY_COLUMNS_A = ["data_", "_pdbx_audit_revision_history.revision_date", "_exptl.method", "_em_3d_reconstruction.resolution", "_refine.ls_d_res_high"]
_NECESSARY_COLUMNS_B = ["data_", "_pdbx_audit_revision_history.revision_date", "_exptl.method", "_em_3d_reconstruction.resolution", "_refine.ls_d_res_high"]
