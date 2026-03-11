"""
Tests for pose clustering (DBSCAN and KMeans).
Uses the synthetic_traj_dir session fixture from conftest.py — no real data files needed.
"""

import pytest
import numpy as np
from sklearn.decomposition import PCA
from pylipid.func import collect_bound_poses, vectorize_poses, write_bound_poses
from pylipid.func import cluster_DBSCAN, cluster_KMeans
from pylipid.api import LipidInteraction
import shutil


@pytest.fixture(scope="module")
def li_and_poses(synthetic_traj_dir, tmp_path_factory):
    save_dir = str(tmp_path_factory.mktemp("cluster_output"))
    li = LipidInteraction(
        synthetic_traj_dir["trajfile_list"],
        cutoffs=synthetic_traj_dir["cutoffs"],
        topfile_list=synthetic_traj_dir["topfile_list"],
        lipid=synthetic_traj_dir["lipid"],
        nprot=1,
        save_dir=save_dir,
    )
    li.collect_residue_contacts()
    li.compute_residue_koff(plot_data=False)
    li.compute_binding_nodes(threshold=2, print_data=False)

    if not li._node_list:
        yield li, {}
    else:
        binding_site_map = {bs_id: nodes for bs_id, nodes in enumerate(li._node_list)}
        pose_pool = collect_bound_poses(
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
        yield li, pose_pool

    shutil.rmtree(save_dir, ignore_errors=True)


class TestCluster:

    def test_cluster_DBSCAN(self, li_and_poses, tmp_path):
        li, pose_pool = li_and_poses
        if not li._node_list:
            pytest.skip("No binding sites detected in synthetic data")
        for bs_id, nodes in enumerate(li._node_list):
            if bs_id not in pose_pool or len(pose_pool[bs_id]) < 5:
                continue
            dist_matrix, pose_traj = vectorize_poses(
                pose_pool[bs_id], nodes, li._protein_ref, li._lipid_ref)
            lipid_dist_per_pose = [dist_matrix[:, pose_id, :].ravel()
                                   for pose_id in np.arange(dist_matrix.shape[1])]
            transformed_data = PCA(n_components=min(2, len(lipid_dist_per_pose) - 1)).fit_transform(lipid_dist_per_pose)
            cluster_labels = cluster_DBSCAN(transformed_data, eps=None, min_samples=None,
                                            metric="euclidean")
            assert len(cluster_labels) == len(lipid_dist_per_pose)
            cluster_id_set = [label for label in np.unique(cluster_labels) if label != -1]
            if cluster_id_set:
                selected_pose_id = [np.random.choice(np.where(cluster_labels == cid)[0], 1)[0]
                                    for cid in cluster_id_set]
                write_bound_poses(pose_traj, selected_pose_id, str(tmp_path),
                                  pose_prefix="BSid{}_cluster_DBSCAN".format(bs_id),
                                  pose_format="gro")

    def test_cluster_KMeans(self, li_and_poses, tmp_path):
        li, pose_pool = li_and_poses
        if not li._node_list:
            pytest.skip("No binding sites detected in synthetic data")
        for bs_id, nodes in enumerate(li._node_list):
            if bs_id not in pose_pool or len(pose_pool[bs_id]) < 5:
                continue
            dist_matrix, pose_traj = vectorize_poses(
                pose_pool[bs_id], nodes, li._protein_ref, li._lipid_ref)
            lipid_dist_per_pose = [dist_matrix[:, pose_id, :].ravel()
                                   for pose_id in np.arange(dist_matrix.shape[1])]
            n_clusters = min(3, len(lipid_dist_per_pose))
            transformed_data = PCA(n_components=min(2, len(lipid_dist_per_pose) - 1)).fit_transform(lipid_dist_per_pose)
            cluster_labels = cluster_KMeans(transformed_data, n_clusters=n_clusters)
            assert len(cluster_labels) == len(lipid_dist_per_pose)
            cluster_id_set = [label for label in np.unique(cluster_labels) if label != -1]
            if cluster_id_set:
                selected_pose_id = [np.random.choice(np.where(cluster_labels == cid)[0], 1)[0]
                                    for cid in cluster_id_set]
                write_bound_poses(pose_traj, selected_pose_id, str(tmp_path),
                                  pose_prefix="BSid{}_cluster_KMeans".format(bs_id),
                                  pose_format="gro")
