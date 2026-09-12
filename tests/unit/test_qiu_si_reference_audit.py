from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from audit_qiu_si_reference import declared_local_imports, safe_member_name


def test_declared_local_imports_only_returns_archived_function_modules():
    source = "from functions_6ref_new2_1 import *\nimport numpy as np\nimport functions_aux\n"
    assert declared_local_imports(source) == ["functions_6ref_new2_1", "functions_aux"]


def test_safe_member_name_rejects_archive_traversal():
    assert safe_member_name("PF_Codes/Polycrystals/functions_2ref.py")
    assert not safe_member_name("../outside.py")
    assert not safe_member_name("/absolute.py")
