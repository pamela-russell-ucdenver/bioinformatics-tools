package tests;

import general.CommandLineParser;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;

import pipeline.ConfigFile;
import pipeline.ConfigFileOption;
import pipeline.ConfigFileSection;


/**
 * @author prussell
 *
 */
public class TestConfigFile {
	
	private static ConfigFileSection sectionCommands = new ConfigFileSection("commands", true);
	private static ConfigFileSection sectionSpecies = new ConfigFileSection("species", false);
	private static ConfigFileSection sectionBasicOptions = new ConfigFileSection("basic_options", true);
	
	private static ConfigFileOption optionSplitTrimBarcodes = new ConfigFileOption("SPLIT_TRIM_BARCODES", 1, false, false, false);
	private static ConfigFileOption optionTrimAdapters = new ConfigFileOption("TRIM_ADAPTERS", 1, false, false, false);
	private static ConfigFileOption optionComputeLibraryStats = new ConfigFileOption("LIBRARY_STATS", 1, false, false, false);
	private static ConfigFileOption optionCountRnaClasses = new ConfigFileOption("RNA_CLASSES", 1, false, false, false);
	private static ConfigFileOption optionFilterRrna = new ConfigFileOption("FILTER_RRNA", 1, false, false, false);
	private static ConfigFileOption optionAlign = new ConfigFileOption("ALIGN", 1, false, false, false);
	private static ConfigFileOption optionAlignToTranscripts = new ConfigFileOption("ALIGN_TO_TRANSCRIPTS", 1, false, false, false);
	private static ConfigFileOption optionRunDge = new ConfigFileOption("RUN_DGE", 1, false, false, false);
	private static ConfigFileOption optionGenomeFasta = new ConfigFileOption("genome_fasta", 2, false, false, true);
	private static ConfigFileOption optionQueueName = new ConfigFileOption("queue_name", 2, false, false, false, "week");
	private static ConfigFileOption optionGenomeBowtieIndex = new ConfigFileOption("genome_bowtie", 2, false, false, true);
	private static ConfigFileOption optionGenomeAssemblyName = new ConfigFileOption("genome_assembly", 2, false, false, true);
	private static ConfigFileOption optionGtfAnnotation = new ConfigFileOption("annotation", 2, false, false, false);
	private static ConfigFileOption optionRnaClassFastaFile = new ConfigFileOption("rna_classes", 3, false, true, false);
	private static ConfigFileOption optionRrnaSequences = new ConfigFileOption("filter_rrna", 2, false, false, false);
	private static ConfigFileOption optionTranscriptFasta = new ConfigFileOption("align_to_transcripts", 2, false, false, false);
	private static ConfigFileOption optionSamplesToMerge = new ConfigFileOption("merge_samples", 10, true, true, false);
	private static ConfigFileOption optionTophatExecutable = new ConfigFileOption("tophat_path", 2, false, false, true);
	private static ConfigFileOption optionBedFileForWig = new ConfigFileOption("bed_file_for_wig", 2, false, false, false);
	private static ConfigFileOption optionWigToBigWigExecutable = new ConfigFileOption("wig_to_bigwig_path", 2, false, false, false);
	private static ConfigFileOption optionWigWriterJar = new ConfigFileOption("wig_writer_jar", 2, false, false, false);
	private static ConfigFileOption optionChrSizeFile = new ConfigFileOption("chr_size_file", 2, false, false, true);
	private static ConfigFileOption optionAlignGlobalStatsJar = new ConfigFileOption("alignment_global_stats_jar_file", 2, false, false, false);
	private static ConfigFileOption optionTophatOption = new ConfigFileOption("tophat_options", 3, true, true, false);
	private static ConfigFileOption optionBowtie2Option = new ConfigFileOption("bowtie2_options", 3, true, true, false);
	private static ConfigFileOption optionNovoalignOption = new ConfigFileOption("novoalign_options", 3, true, true, false);
	private static ConfigFileOption optionFragmentSizeDistOption = new ConfigFileOption("fragment_size_dist_options", 3, true, true, false);
	private static ConfigFileOption optionBowtie2Executable = new ConfigFileOption("bowtie2_path", 2, false, false, true);
	private static ConfigFileOption optionBowtie2BuildExecutable = new ConfigFileOption("bowtie2_build_path", 2, false, false, true);
	private static ConfigFileOption optionSamtoolsPath = new ConfigFileOption("samtools_path", 2, false, false, true);
	private static ConfigFileOption optionIgvToolsExecutable = new ConfigFileOption("igvtools_path", 2, false, false, true);
	private static ConfigFileOption optionFastqReadNumberDelimiter = new ConfigFileOption("fastq_read_number_delimiter", 2, false, false, false, "\\s++");
	private static ConfigFileOption optionPicardDir = new ConfigFileOption("picard_directory", 2, false, false, true);
	private static ConfigFileOption optionFastxDirectory = new ConfigFileOption("fastx_directory", 2, false, false, false);
	private static ConfigFileOption optionRead1Adapter = new ConfigFileOption("sequencing_adapter_read1", 2, false, false, false);
	private static ConfigFileOption optionRead2Adapter = new ConfigFileOption("sequencing_adapter_read2", 2, false, false, false);
	private static ConfigFileOption optionNovoalignExecutable = new ConfigFileOption("novoalign_path", 2, false, false, false);
	private static ConfigFileOption optionGenomeNovoindex = new ConfigFileOption("genome_novoindex", 2, false, false, false);
	private static ConfigFileOption optionPicardRefFlat = new ConfigFileOption("picard_ref_flat", 2, false, false, false);
	private static ConfigFileOption optionPicardRibosomalIntervals = new ConfigFileOption("picard_ribosomal_intervals", 2, false, false, false);
	private static ConfigFileOption optionPicardStrandSpecificity = new ConfigFileOption("picard_strand_specificity", 2, false, false, false);
	private static ConfigFileOption optionPicardCollectRnaSeqMetricsJar = new ConfigFileOption("picard_collect_rnaseq_metrics", 2, false, false, false);
	private static ConfigFileOption optionTranscriptomeSpaceStatsBedFile = new ConfigFileOption("bed_transcriptome_space_stats", 2, false, false, false);
	private static ConfigFileOption optionGenomicSpaceStatsSizeFile = new ConfigFileOption("size_file_genomic_space_stats", 2, false, false, false);
	private static ConfigFileOption optionSpeciesName = new ConfigFileOption("species", 2, false, false, false, "mouse");

	
	private static ConfigFile getConfigFile(String fileName) throws IOException {
		
		sectionCommands.addAllowableOption(optionAlign);
		sectionCommands.addAllowableOption(optionAlignToTranscripts);
		sectionCommands.addAllowableOption(optionComputeLibraryStats);
		sectionCommands.addAllowableOption(optionCountRnaClasses);
		sectionCommands.addAllowableOption(optionFilterRrna);
		sectionCommands.addAllowableOption(optionRunDge);
		sectionCommands.addAllowableOption(optionSplitTrimBarcodes);
		sectionCommands.addAllowableOption(optionTrimAdapters);
		
		sectionSpecies.addAllowableOption(optionSpeciesName);
		
		sectionBasicOptions.addAllowableOption(optionGenomeFasta);
		sectionBasicOptions.addAllowableOption(optionQueueName);
		sectionBasicOptions.addAllowableOption(optionGenomeBowtieIndex);
		sectionBasicOptions.addAllowableOption(optionGenomeAssemblyName);
		sectionBasicOptions.addAllowableOption(optionGtfAnnotation);
		sectionBasicOptions.addAllowableOption(optionRnaClassFastaFile);
		sectionBasicOptions.addAllowableOption(optionRrnaSequences);
		sectionBasicOptions.addAllowableOption(optionTranscriptFasta);
		sectionBasicOptions.addAllowableOption(optionSamplesToMerge);
		sectionBasicOptions.addAllowableOption(optionTophatExecutable);
		sectionBasicOptions.addAllowableOption(optionBedFileForWig);
		sectionBasicOptions.addAllowableOption(optionWigToBigWigExecutable);
		sectionBasicOptions.addAllowableOption(optionWigWriterJar);
		sectionBasicOptions.addAllowableOption(optionChrSizeFile);
		sectionBasicOptions.addAllowableOption(optionAlignGlobalStatsJar);
		sectionBasicOptions.addAllowableOption(optionTophatOption);
		sectionBasicOptions.addAllowableOption(optionBowtie2Option);
		sectionBasicOptions.addAllowableOption(optionNovoalignOption);
		sectionBasicOptions.addAllowableOption(optionFragmentSizeDistOption);
		sectionBasicOptions.addAllowableOption(optionBowtie2Executable);
		sectionBasicOptions.addAllowableOption(optionBowtie2BuildExecutable);
		sectionBasicOptions.addAllowableOption(optionSamtoolsPath);
		sectionBasicOptions.addAllowableOption(optionIgvToolsExecutable);
		sectionBasicOptions.addAllowableOption(optionFastqReadNumberDelimiter);
		sectionBasicOptions.addAllowableOption(optionPicardDir);
		sectionBasicOptions.addAllowableOption(optionFastxDirectory);
		sectionBasicOptions.addAllowableOption(optionRead1Adapter);
		sectionBasicOptions.addAllowableOption(optionRead2Adapter);
		sectionBasicOptions.addAllowableOption(optionNovoalignExecutable);
		sectionBasicOptions.addAllowableOption(optionGenomeNovoindex);
		sectionBasicOptions.addAllowableOption(optionPicardRefFlat);
		sectionBasicOptions.addAllowableOption(optionPicardRibosomalIntervals);
		sectionBasicOptions.addAllowableOption(optionPicardStrandSpecificity);
		sectionBasicOptions.addAllowableOption(optionPicardCollectRnaSeqMetricsJar);
		sectionBasicOptions.addAllowableOption(optionTranscriptomeSpaceStatsBedFile);
		sectionBasicOptions.addAllowableOption(optionGenomicSpaceStatsSizeFile);
		
		Collection<ConfigFileSection> sections = new ArrayList<ConfigFileSection>();
		sections.add(sectionCommands);
		sections.add(sectionSpecies);
		sections.add(sectionBasicOptions);
		
		ConfigFile cf = new ConfigFile(sections, fileName);
		return cf;

	}

	
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
		CommandLineParser p = new CommandLineParser();
		p.addStringArg("-c", "Config file", true);
		p.parse(args);
		String configFile = p.getStringArg("-c");
		
		@SuppressWarnings("unused")
		ConfigFile cf = getConfigFile(configFile);
		
		
		
	}

}
