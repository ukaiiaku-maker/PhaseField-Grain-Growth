from pathlib import Path
import numpy as np,pytest
from scripts.prepare_qiu_control_segment import CONTROLS,step
def test_control_codes_are_frozen():assert CONTROLS=={"QIU_SI_4REF_A0_PORT":0,"QIU_SI_4REF_ANISO_E":1,"QIU_SI_4REF_ANISO_M":2,"QIU_SI_4REF_ANISO_EM_INV":3}
def test_reads_exact_parent_step(tmp_path):
 p=tmp_path/'p.npz';np.savez(p,accepted_step=np.asarray(1000));assert step(p)==1000
