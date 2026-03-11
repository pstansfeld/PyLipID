"""
Tests for binding site detection, pose collection, scoring and surface area.
Uses the synthetic_traj_dir session fixture from conftest.py — no real data files needed.
"""

import pytest
import numpy as np
from pylipid.func import (get_node_list, collect_bound_poses, vectorize_poses,
                           calculate_scores, write_bound_poses)
from pylipid.api import LipidInteraction
import shutil


@pytest.fixture(scope="module")
def li(synthetic_traj_dir, tmp_path_factory):
    save_dir = str(tmp_path_factory.mktemp("bs_output"))
    instance = LipidInteraction(
        synthetic_traj_dir["trajfile_list"],
        cutoffs=synthetic_traj_dir["cutoffs"],
        topfile_list=synthetic_traj_dir["topfile_list"],
        lipid=synthetic_traj_dir["lipid"],
        nprot=1,
        save_dir=save_dir,
    )
    instance.collect_residue_contacts()
    instance.compute_residue_koff(plot_data=False)
    instance.compute_binding_nodes(threshold=2, print_data=False)
    yield instance
    shutil.rmtree(save_dir, ignore_errors=True)


@pytest.fixture(scope="module")
def pose_pool(li, synthetic_traj_dir):
    """Collect bound poses for all detected binding sites."""
    if not li._node_list:
        pytest.skip("No binding sites detected in synthetic data")
    binding_site_map = {bs_id: nodes for bs_id, nodes in enumerate(li._node_list)}
    return collect_bound_poses(
        binding_site_map,
        li._contact_residues_low,
        li.trajfile_list,
        li.topfile_list,
        li.lipid,
        li._protein_ref,
        li._lipid_ref,
        stride=li.stride,
        nprot=li.nprot,
    )


class TestBindingSites:

    def test_get_node_list(self, li):
        corrcoef = li.interaction_corrcoef
        node_list, modularity = get_node_list(corrcoef, threshold=2)
        assert isinstance(node_list, list)

    def test_collect_binding_poses(self, pose_pool):
        assert isinstance(pose_pool, dict)


class TestBindingPoses:

    def test_vectorize_poses(self, li, pose_pool):
        for bs_id, nodes in enumerate(li._node_list):
            if bs_id not in pose_pool or len(pose_pool[bs_id]) == 0:
                continue
            dist_matrix, pose_traj = vectorize_poses(
                pose_pool[bs_id], nodes, li._protein_ref, li._lipid_ref)
            assert dist_matrix.shape[0] == li._lipid_ref.n_atoms
            assert dist_matrix.shape[1] == len(pose_pool[bs_id])
            assert dist_matrix.shape[2] == len(nodes)

    def test_calculate_scores(self, li, pose_pool):
        for bs_id, nodes in enumerate(li._node_list):
            if bs_id not in pose_pool or len(pose_pool[bs_id]) == 0:
                continue
            dist_matrix, pose_traj = vectorize_poses(
                pose_pool[bs_id], nodes, li._protein_ref, li._lipid_ref)
            scores = calculate_scores(dist_matrix)
            assert len(scores) == pose_traj.n_frames
            scores_weighted = calculate_scores(dist_matrix, score_weights={"RHO": 10})
            assert len(scores_weighted) == pose_traj.n_frames

    def test_write_binding_poses(self, li, pose_pool, tmp_path):
        for bs_id, nodes in enumerate(li._node_list):
            if bs_id not in pose_pool or len(pose_pool[bs_id]) == 0:
                continue
            dist_matrix, pose_traj = vectorize_poses(
                pose_pool[bs_id], nodes, li._protein_ref, li._lipid_ref)
            scores = calculate_scores(dist_matrix)
            if len(scores) == 0:
                continue
            num_of_poses = min(3, pose_traj.n_frames)
            pose_indices = np.argsort(scores)[::-1][:num_of_poses]
            write_bound_poses(pose_traj, pose_indices, str(tmp_path),
                              pose_prefix="BSid{}_top".format(bs_id),
                              pose_format="gro")

    def test_compute_surface_area(self, li):
        # Use the high-level API method — calculate_surface_area_wrapper is internal
        li.compute_surface_area()
