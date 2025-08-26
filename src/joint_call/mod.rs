mod find_inversions;
mod get_refined_svs;
mod joint_call_all_samples;
mod merge_cnv;
mod merge_haplotypes;
mod read_sample_data;
mod sample_output;
mod supporting_read_names;
mod write_merged_contigs;

use self::find_inversions::find_inversions;
use self::joint_call_all_samples::joint_genotype_all_samples;
use self::merge_cnv::merge_sv_with_depth_info;
use self::merge_haplotypes::merge_haplotypes;
pub use self::read_sample_data::SampleJointCallData;
use self::read_sample_data::{SharedJointCallData, read_all_sample_data};
use self::sample_output::setup_sample_output;

use crate::cli::{JointCallDerivedSettings, JointCallSettings, SharedSettings};
use crate::globals::PROGRAM_VERSION;
use crate::large_variant_output::{VcfSettings, write_indexed_sv_vcf_file};
use crate::run_stats::{JointCallRunStats, ScoreStats, RunStep, MultiSampleMergeStats, delete_run_stats, write_joint_call_run_stats};
use crate::sv_group::{SVGroup};


use log::info;
//use std::fs::File;
use serde::{Serialize, de::DeserializeOwned};
use unwrap::unwrap;


/// Serialize data to MessagePack format and write to file
pub fn serialize_to_mpack<T: Serialize>(filename: &str, data: &T) {
    let mut buf = Vec::new();
    data.serialize(&mut rmp_serde::Serializer::new(&mut buf))
         .unwrap();

    info!("Writing data to binary file: '{filename}'");
    unwrap!(
            std::fs::write(&filename, buf.as_slice()),
            "Unable to open and write binary file: '{filename}'"
        );
}

/// Deserialize data from MessagePack format file
pub fn deserialize_from_mpack<T: DeserializeOwned>(filename: &str) -> T {
    let buf = unwrap!(
            std::fs::read(&filename),
            "Unable to open SV groups binary file: '{filename}'"
        );
        unwrap!(
            rmp_serde::from_slice(&buf),
            "Unable to parse SV groups binary file: '{filename}'"
        )
}


fn merge_shards(settings: &JointCallSettings) -> (Vec<SVGroup>, ScoreStats) {
    let mut all_vectors: Vec<Vec<SVGroup>> = Vec::new();
    let mut combined_stats = ScoreStats::default();
    for chunk in 0..=settings.chunk_size { //TODO: handle chunk numbers

        let content =deserialize_from_mpack::<Vec<SVGroup>>(&[&settings.output_dir.to_string(),"SV_genotyped_",&chunk.to_string(),".bin"].join(""));
        all_vectors.push(content);
        let stats =deserialize_from_mpack::<ScoreStats>(&[&settings.output_dir.to_string(),"SV_genotyped_stats_",&chunk.to_string(),".bin"].join(""));
        combined_stats.merge(&stats);
    }

    // Flatten into a single vector
   (all_vectors.into_iter().flatten().collect(), combined_stats)

}

pub fn run_joint_call(
    shared_settings: &SharedSettings,
    settings: &JointCallSettings,
    derived_settings: &JointCallDerivedSettings,
) {
    // Now that we're committed to a run, remove any possible older run stats file that could be present in case this is a clobber run
    //
    // The run stats file is used as a marker of a successfully finished run, so removing it here allows run completion to be determined
    // from whether the new file is written at the end of this discover step.
    //
    delete_run_stats(&settings.output_dir);

    let (shared_data, all_sample_data) =
        read_all_sample_data(shared_settings, settings, derived_settings);

    if settings.stages.contains(&"bigwig".to_string()) {
        setup_sample_output(
            &settings.output_dir,
            &shared_data.chrom_list,
            &shared_data,
            &all_sample_data,
        );
    }
    let mut merged_sv_groups; //: Vec<SVGroup> = Vec::new();
    let mut merge_stats; // = MultiSampleMergeStats::default();

    if settings.stages.contains(&"merge".to_string()) {
        (merged_sv_groups, merge_stats) =
            merge_haplotypes(shared_settings, settings, &shared_data, &all_sample_data);
        
        write_merged_contigs::write_merged_contig_alignments(
            shared_settings,
            settings,
            &shared_data,
            &merged_sv_groups,
        );

        serialize_to_mpack(&[&settings.output_dir.to_string(),"/SV_clusters.bin"].join(""), &merged_sv_groups);
        serialize_to_mpack(&[&settings.output_dir.to_string(),"/SV_clusters_stats.bin"].join(""), &merge_stats);

    } else {
        merged_sv_groups = deserialize_from_mpack::<Vec<SVGroup>>(&[&settings.output_dir.to_string(),"/SV_clusters.bin"].join(""));
        merge_stats = deserialize_from_mpack::<MultiSampleMergeStats>(&[&settings.output_dir.to_string(),"/SV_clusters_stats.bin"].join(""));
    }

    //skip any chunking if chunk_size=1
    if settings.chunk_size > 1 {
        let total_size: usize = merged_sv_groups.len();
        println!("Chunking SV groups (n={}) into {} parts", total_size, total_size / settings.chunk_size);
        //TODO: how do we handle perfect file spliting vs remainder
       //merged_sv_groups = merged_sv_groups.chunks(total_size / settings.chunk_size).nth(settings.chunk).unwrap().to_vec();

        if merged_sv_groups.len() % settings.chunk_size > 0 && settings.chunk == (merged_sv_groups.len() / settings.chunk_size) { //handle the remainder
            merged_sv_groups = merged_sv_groups.chunks_exact(total_size / settings.chunk_size).remainder().to_vec();
        } else {
            merged_sv_groups = merged_sv_groups.chunks_exact(total_size / settings.chunk_size).nth(settings.chunk).unwrap().to_vec();
        }
    }

    let enable_phasing = true;
    let (mut sv_groups, mut score_stats) = joint_genotype_all_samples(
        shared_settings,
        settings,
        enable_phasing,
        &shared_data,
        &all_sample_data,
        merged_sv_groups,
    );
    
    if settings.chunk_size > 1  && settings.stages.contains(&"genotype".to_string()) {
    serialize_to_mpack(&[&settings.output_dir.to_string(),"SV_genotyped_",&settings.chunk.to_string(),".bin"].join(""), &sv_groups);
    serialize_to_mpack(&[&settings.output_dir.to_string(),"SV_genotyped_stats_",&settings.chunk.to_string(),".bin"].join(""), &score_stats);
    }

    //merge shared vectors
    //only do this if chunk_size > 1
    
    if settings.stages.contains(&"finish".to_string()) {
        if settings.chunk_size > 1  {
            (sv_groups, score_stats) = merge_shards(settings);
        }

    find_inversions(&mut sv_groups);

    let refined_cnvs =
        merge_sv_with_depth_info(settings, &shared_data, &all_sample_data, &mut sv_groups);

    let vcf_settings = VcfSettings::new(
        &shared_data.ref_filename,
        &settings.output_dir,
        settings.min_qual,
        settings.no_vcf_dedup,
        enable_phasing,
        settings.treat_single_copy_as_haploid,
    );

    let sample_names = all_sample_data
        .iter()
        .map(|x| x.sample_name.as_str())
        .collect::<Vec<_>>();
    let vcf_stats = write_indexed_sv_vcf_file(
        shared_settings,
        &vcf_settings,
        &shared_data.genome_ref,
        &shared_data.chrom_list,
        &sample_names,
        &sv_groups,
        &refined_cnvs,
        false,
    );

    score_stats.vcf_output_record_count = vcf_stats.output_record_count;
    score_stats.vcf_duplicate_record_count = vcf_stats.duplicate_record_count;

    if settings.report_supporting_reads {
        supporting_read_names::write_supporting_read_names(
            &settings.output_dir,
            &sample_names,
            &sv_groups,
        );
    }
    // can remove all the mpack files
    // In addition to useful statistics this file acts as a marker for a successfully completed run, so it must be written last.
    write_joint_call_run_stats(
        &settings.output_dir,
        &JointCallRunStats {
            run_step: RunStep {
                name: "joint-call".to_string(),
                version: PROGRAM_VERSION.to_string(),
            },
            merge_stats,
            score_stats,
        },
    );
    }
}
