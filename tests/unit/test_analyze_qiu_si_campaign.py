import io
import json
import tarfile

import numpy as np

from scripts.analyze_qiu_si_campaign import analyze


def make_archive(path, step):
    phi=np.zeros((2,4,4));phi[0,:,:2]=1;phi[1,:,2:]=1
    labels=np.argmax(phi,axis=0);delta=np.zeros_like(phi);delta[0,0,1]=.01;delta[1,0,1]=-.01
    dtype=[('line_density_x','f8'),('line_density_y','f8')];lines=np.zeros(2,dtype=dtype);lines['line_density_x']=[1,-1]
    summary={"grain_count":2,"boundary_density":.5,"compactness_mean":.7,"active_support_mean":1.0,"active_support_max":1,"phase_sum_max_abs_error":0.0}
    buffer=io.BytesIO()
    np.savez_compressed(buffer,accepted_step=np.asarray(step),physical_time=np.asarray(step*.1),output_scaling=np.asarray([1.,1.]),phi=phi,labels=labels,orientations=np.asarray([0.,.2]),accepted_delta=delta,sigma11=np.zeros((4,4)),sigma12=np.ones((4,4)),sigma22=np.zeros((4,4)),gb_records=lines,pair_cell_x=np.asarray([0]),pair_cell_y=np.asarray([1]),pair_phase_i=np.asarray([0]),pair_phase_j=np.asarray([1]),signed_elastic_pair_force=np.asarray([2.]),capillary_pair_term=np.asarray([3.]),barrier_pair_term=np.asarray([4.]),pair_beta_first=np.asarray([.1]),pair_reference_first=np.asarray([0.]),pair_reference_second=np.asarray([.785]),summary_json=np.asarray(json.dumps(summary)))
    payload=buffer.getvalue();info=tarfile.TarInfo(f"./output/checkpoints/native-state-step{step:06d}.npz");info.size=len(payload)
    with tarfile.open(path,"w:gz") as stream:stream.addfile(info,io.BytesIO(payload))


def test_analysis_builds_machine_readable_tables_figures_and_movies(tmp_path):
    native=tmp_path/"native.tar.gz";aniso=tmp_path/"aniso.tar.gz";make_archive(native,250);make_archive(aniso,250)
    manifest=tmp_path/"manifest.json";manifest.write_text(json.dumps({"cases":[{"case":"QIU_SI_NATIVE","archives":[str(native)],"movie":True},{"case":"QIU_SI_ANISO_EM_INV","archives":[str(aniso)],"movie":True}]}))
    output=tmp_path/"analysis";summary=analyze(manifest,output)
    assert summary["classification"]=="QIU_SI_FINAL_ANALYSIS_COMPLETE"
    assert summary["matched_native_anisotropic_states"]==1
    assert (output/"matched_time_native_anisotropic.csv").is_file()
    assert (output/"figures/energy_work_balance.pdf").is_file()
    assert (output/"movies/qiu_si_aniso_em_inv.gif").is_file()
    assert (output/"checksums.sha256").is_file()
