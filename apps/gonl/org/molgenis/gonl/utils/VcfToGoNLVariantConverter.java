package org.molgenis.gonl.utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.codec.digest.DigestUtils;
import org.apache.log4j.BasicConfigurator;
import org.apache.log4j.Logger;
import org.molgenis.core.OntologyTerm;
import org.molgenis.pheno.ObservableFeature;
import org.molgenis.pheno.ObservedValue;
import org.molgenis.pheno.Panel;
import org.molgenis.pheno.Species;
import org.molgenis.util.CsvFileWriter;
import org.molgenis.util.CsvWriter;
import org.molgenis.util.Entity;
import org.molgenis.util.vcf.VcfReader;
import org.molgenis.util.vcf.VcfReaderListener;
import org.molgenis.util.vcf.VcfRecord;
import org.molgenis.variant.Chromosome;
import org.molgenis.variant.GenomeBuild;
import org.molgenis.variant.SequenceVariant;

/**
 * This method can extract cohort level or individual level data from vcf files
 * and convert this into a pheno model compatible data set.
 * 
 * Currently, only aggregate data is exported as observedvalue.
 */
public class VcfToGoNLVariantConverter
{
	private static final Logger logger = Logger.getLogger(VcfToGoNLVariantConverter.class);

	public static final int BATCH_SIZE = 100000;

	public static void main(String[] args) throws Exception
	{

		if (args.length != 4)
		{
			// /Users/despoina/Documents/_____VCFpathoData/vcf_files/gonl.chr22.release4.sum.vcf
			// /Users/despoina/Documents/_____VCFpathoData/output/ test test
			System.err.println("Usage: <input.vcf> <output.dir> <panel_name> <encryption_salt>");
			return;
		}

		BasicConfigurator.configure();

		File vcfFile = new File(args[0]);
		if (!vcfFile.exists()) throw new FileNotFoundException("input does not exist: " + args[0]);
		else if (vcfFile.isDirectory()) throw new IOException("input is a directory: " + args[0]);
		File outputDir = new File(args[1]);
		if (!outputDir.exists()) if (!outputDir.mkdir()) throw new IOException("Could not create directory: "
				+ outputDir);
		else if (!outputDir.isDirectory()) throw new IOException("output directory is not a directory");

		VcfToGoNLVariantConverter convert = new VcfToGoNLVariantConverter(args[2], args[3]);
		convert.convertVariants(vcfFile, outputDir);
	}

	private final String panelName;
	/**
	 * Random data used as additional input for individual name encryption
	 */
	private final String salt;

	public VcfToGoNLVariantConverter(String panelName, String salt)
	{
		if (panelName == null) throw new IllegalArgumentException();
		if (salt == null) throw new IllegalArgumentException();
		this.panelName = panelName;
		this.salt = salt;
	}

	public void convertVariants(final File vcfFile, final File outputDir) throws Exception
	{
		System.out.println("converting aggregate data from vcf=" + vcfFile + " to directory " + outputDir);
		final VcfReader vcf = new VcfReader(vcfFile);

		SimpleDateFormat sdf = new SimpleDateFormat("MMMM d, yyyy");
		String date = sdf.format(new Date());

		// final List<Variant> variants = new ArrayList<Variant>();
		final List<SequenceVariant> variants = new ArrayList<SequenceVariant>();
		final List<ObservedValue> values = new ArrayList<ObservedValue>();
		final List<String> chromosomes = new ArrayList<String>();
		final List<String> dbXrefs = new ArrayList<String>();
		final List<GenomeBuild> builds = new ArrayList<GenomeBuild>();
		final List<Species> species = new ArrayList<Species>();
		final List<ObservableFeature> features = new ArrayList<ObservableFeature>();
		final Panel panel = new Panel();

		// create file names
		final File fileVariants = new File(outputDir.getAbsolutePath() + File.separatorChar
				+ SequenceVariant.class.getSimpleName() + ".txt");
		final File fileObservedValues = new File(outputDir.getAbsolutePath() + File.separatorChar
				+ ObservedValue.class.getSimpleName() + ".txt");

		// create file headers
		final String[] variantHeaders = new String[]
		{ SequenceVariant.ID, SequenceVariant.NAME, SequenceVariant.DESCRIPTION, SequenceVariant.INVESTIGATION_NAME,
				SequenceVariant.ONTOLOGYREFERENCE_NAME, SequenceVariant.ALTERNATEID_NAME, SequenceVariant.__TYPE,
				SequenceVariant.ALT, SequenceVariant.CHR, SequenceVariant.CHR_NAME, SequenceVariant.DBREFS,
				SequenceVariant.DBREFS_NAME, SequenceVariant.ENDBP, SequenceVariant.ENDBP, SequenceVariant.REF,
				SequenceVariant.STARTBP, SequenceVariant.VARIANTTYPE };

		final String[] ovHeaders = new String[]
		{ ObservedValue.ID, ObservedValue.INVESTIGATION_NAME, ObservedValue.PROTOCOLAPPLICATION_NAME,
				ObservedValue.FEATURE_NAME, ObservedValue.TARGET_NAME, ObservedValue.ONTOLOGYREFERENCE_NAME,
				ObservedValue.VALUE, ObservedValue.CHARACTERISTICVALUE_NAME, ObservedValue.CHARACTERISTICVALUE_NAME,
				ObservedValue.RELATION_NAME, ObservedValue.TIME, ObservedValue.ENDTIME };
		// create files
		// createFileAndHeader(fileVariants, variantHeaders);
		// createFileAndHeader(fileObservedValues, ovHeaders);

		final List<Integer> count = new ArrayList<Integer>();
		count.add(0);

		final String encryptionSalt = this.salt;
		final Map<String, String> encIndividualMap = new HashMap<String, String>();
		vcf.parse(new VcfReaderListener()
		{
			@Override
			public void handleLine(int lineNumber, VcfRecord record) throws Exception
			{
				// create hgvs notation, one record per variant reported (even
				// if they
				// are on same location)
				List<String> alt = record.getAlt();
				for (int i = 0; i < alt.size(); i++)
				{
					SequenceVariant v = new SequenceVariant();
					ObservedValue observedValue = new ObservedValue();
					List<ObservedValue> observedValuesList = new ArrayList<ObservedValue>();

					// TODO result is not used, do we need this?
					// String result = "chr" + record.getChrom() + ":g.";

//#CHROM 	POS		ID			REF	ALT			QUAL	FILTER									INFO	
//22	 	12345	.			A	T			2000	TruthSensitivityTranche99.60to99.70		AC=3;AN=10;GTC=3,1,1	
//22		12346	rs0987654	A	G,C			2000	TruthSensitivityTranche99.60to99.70		AC=1,4;AN=10;GTC=2,0,0,1,1,1	
//22		12347	rs0987654	A	G,C,TGCCAAT	2000	TruthSensitivityTranche99.60to99.70		AC=1,6,3;AN=20;GTC=4,0,0,2,1,1,0,0,1,1

					v.setName("chr" + record.getChrom() + ":g." + record.getPos() + record.getRef() + ">" + alt.get(i));
					v.setName(v.getName().replace("|", "_"));

					// ref
					v.setRef(record.getRef());

					// alt
					v.setAlt(alt.get(i));

					// ref TODO : check this is used in Variant, do we need to
					// set something in SequeceVariant?
					// v.setResidues(record.getRef());

					// alt TODO : check this is used in Variant, do we need to
					// set something in SequeceVariant?
					// v.setAltResidues(alt.get(i));

					// chr
					String chromName = record.getChrom().replace("|", "_");
					v.setChr_Name(chromName);

					// check if chrom exists, otherwise add
					System.out.println(">>>>" + chromosomes);
					if (!chromosomes.contains(chromName)) chromosomes.add(chromName);

					// pos
					v.setStartBP(record.getPos());
					v.setEndBP(record.getPos());

					// pos TODO : pick one of these two, Pieter?
					// TODO : check this is used in Variant, do we need to set
					// something in SequeceVariant?
					// v.setStartGdna(record.getPos());
					// v.setEndGdna(record.getPos());

					v.setDescription("" + record.getId());

					// put alt allele counts in description
					System.out.println(vcf.getInfoFields());
					for (int j = 0; j != vcf.getInfoFields().size(); j++)
					{
						System.out.println(vcf.getInfoFields().get(j));
						if (record.getInfo(vcf.getInfoFields().get(j)) != null)
						{
							// System.out.println(record.getInfo(vcf.getInfoFields().get(j)).get(i));
							List<String> var3 = vcf.getInfoFields();
							String key = var3.get(j);
							List<String> info = record.getInfo(key);
							if (info == null | info.isEmpty()) observedValue.setValue(info.get(0));

							else
								logger.warn("unknown key: " + key);

						}

//##INFO=<ID=GTC,Number=G,Type=Integer,Description="GenoType Counts. For each ALT allele in the same order as listed = 0/0,0/1,1/1,0/2,1/2,2/2,0/3,1/3,2/3,3/3,etc. 
//Phasing is ignored; hence 1/0, 0|1 and 1|0 are all counted as 0/1. When one or more alleles is not called for a genotype in a specific sample (./., ./0, ./1, ./2, etc.),
//that sample's genotype is completely discarded for calculating GTC.">

						
						// create 3 features : AC, AN, GTC
						ObservableFeature acf = new ObservableFeature();
						acf.setName("AC");
						features.add(acf);

						ObservableFeature anf = new ObservableFeature();
						acf.setName("AN");
						features.add(anf);

						ObservableFeature gtcf = new ObservableFeature();
						acf.setName("GTC");
						features.add(gtcf);

						ObservedValue ac = new ObservedValue();
						ac.setFeature(acf);
						ac.setFeature_Name("AC");
						ac.setValue(record.getInfo("AC").get(i));

						ObservedValue an = new ObservedValue();
						ac.setFeature(anf);
						an.setFeature_Name("AN");
						an.setValue(record.getInfo("AN").get(i));

						ObservedValue gtc = new ObservedValue();
						ac.setFeature(gtcf);
						gtc.setFeature_Name("GTC");
						// gtc.setValue(record.getInfo("GTC").get(i));

						observedValuesList.add(ac);
						observedValuesList.add(an);
						observedValuesList.add(gtc);

						observedValue.setValue(record.getInfo("AC").get(i));
						observedValue.setValue(record.getInfo("AN").get(i));
						// GTC is not only one number
						String tmp = "";
						for (int k = 0; k != record.getInfo("GTC").size(); k++)
							tmp = tmp + record.getInfo("GTC").get(k) + " , ";
						observedValue.setValue(tmp);

					}

					variants.add(v);
					values.add(observedValue);

					observedValue.setRelation_Name("Allele count");
					observedValue.setFeature_Name(v.getName());

				}
				List<ObservedValue> observedValuesList = new ArrayList<ObservedValue>();

				for (ObservedValue ov : observedValuesList)
				{
					values.add(ov);
				}

				if (variants.size() >= BATCH_SIZE)
				{
					writeBatch(variants, fileVariants, variantHeaders);
					variants.clear();
					writeBatch(values, fileObservedValues, ovHeaders);
					values.clear();

					count.set(0, count.get(0) + BATCH_SIZE);

					System.out.println(new Date() + " converted variants:" + count.get(0));
				}
			}
		});

		// write remaining data for last batch.
		writeBatch(variants, fileVariants, variantHeaders);
		writeBatch(values, fileObservedValues, ovHeaders);

		// write chromsomes
		List<Chromosome> chrList = new ArrayList<Chromosome>();
		int order = 0;
		for (String chr : chromosomes)
		{
			Chromosome c = new Chromosome();
			c.setName(chr);
			c.setGenomeBuild_Name("hg19");
			c.setOrderNr(++order);
			// this was in GONL before merge
			// if (chr.matches("[a-zA-Z]+"))
			// {
			// c.setIsAutosomal(false);
			// }
			if (chr.matches("[a-zA-Z]+"))
			{
				c.setIsAutosomal("no");
			}
			// else if (chr.matches("^CH?R?"))
			// {
			// c.setIsAutosomal("no");
			// }
			else if (chr.startsWith("UN"))
			{
				c.setIsAutosomal("unknown");
			}
			else
			{
				c.setIsAutosomal("yes");
			}

			chrList.add(c);
		}

		// write dbXrefs
		List<OntologyTerm> ontoList = new ArrayList<OntologyTerm>();

		// add ontlogy terms
		Species s = new Species();
		s.setName("homo sapiens");
		species.add(s);

		File chrFile = new File(outputDir.getAbsolutePath() + File.separatorChar + Chromosome.class.getSimpleName()
				+ ".txt");
		String[] chrHeader = new String[]
		// in GONL before merge { Chromosome.NAME, Chromosome.GENOMEBUILD_NAME,
		// Chromosome.ORDERNR, Chromosome.ISAUTOSOMAL };
		{ Chromosome.ID, Chromosome.NAME, Chromosome.DESCRIPTION, Chromosome.INVESTIGATION_NAME,
				Chromosome.ONTOLOGYREFERENCE_NAME, Chromosome.ALTERNATEID_NAME, Chromosome.BPLENGTH,
				Chromosome.GENOMEBUILD, Chromosome.ORDERNR, Chromosome.ISAUTOSOMAL, Chromosome.BPLENGTH,
				Chromosome.GENOMEBUILD_NAME };
		// {
		createFileAndHeader(chrFile, chrHeader);
		writeBatch(chrList, chrFile, chrHeader);

		File ontoFile = new File(outputDir.getAbsolutePath() + File.separatorChar + OntologyTerm.class.getSimpleName()
				+ ".txt");
		String[] ontoHeader = new String[]
		// GONL before merger { "name" };
		{ "id", "name", "ontology_name", "termAccession", "definition", "termPath" };
		createFileAndHeader(ontoFile, ontoHeader);
		writeBatch(ontoList, ontoFile, ontoHeader);

		GenomeBuild build = new GenomeBuild();
		build.setSpecies_Name("homo sapiens");
		build.setName("hg19");
		builds.add(build);

		final File buildsFile = new File(outputDir.getAbsolutePath() + File.separatorChar
				+ GenomeBuild.class.getSimpleName() + ".txt");
		// final File buildsFile = new File(targetDir, "GenomeBuild.txt");

		String[] buildsHeader = new String[]
		{ "name", "species_name" };
		createFileAndHeader(buildsFile, buildsHeader);
		writeBatch(builds, buildsFile, buildsHeader);

		panel.setName(this.panelName);

		// final File panelFile = new File(outputDir.getAbsolutePath() +
		// File.separatorChar + Panel.class.getSimpleName()
		// + ".txt");
		// String[] panelHeader = new String[]
		// { "name", "individuals_name" };
		// createFileAndHeader(panelFile, panelHeader);
		// writeBatch(Collections.singletonList(panel), panelFile, panelHeader);

		ObservableFeature f = new ObservableFeature();
		f.setName("Allele count");
		features.add(f);

		final File featureFile = new File(outputDir.getAbsolutePath() + File.separatorChar
				+ ObservableFeature.class.getSimpleName() + ".txt");
		String[] featureHeader = new String[]
		{ "name" };
		createFileAndHeader(featureFile, featureHeader);
		writeBatch(features, featureFile, featureHeader);

		final File speciesFile = new File(outputDir.getAbsolutePath() + File.separatorChar
				+ Species.class.getSimpleName() + ".txt");
		String[] speciesHeader = new String[]
		{ "id", "name", "ontology_name", "termAccession", "definition", "termPath" };
		// final File speciesFile = new File(targetDir, "Species.txt");
		createFileAndHeader(speciesFile, speciesHeader);
		writeBatch(species, speciesFile, speciesHeader);

	}

	private String encrypt(String value, String salt)
	{
		return DigestUtils.md5Hex(value + salt);
	}

	private void createFileAndHeader(File file, String[] fields) throws IOException
	{
		CsvWriter writer = new CsvFileWriter(file, Arrays.asList(fields));
		writer.setSeparator(',');
		try
		{
			writer.writeHeader();
		}
		finally
		{
			writer.close();
		}
	}

	private void writeBatch(List<? extends Entity> entities, File file, String[] fields) throws IOException
	{
		if (entities.size() > 0)
		{
			System.out.println("Writing to " + file);

			// create appending csvWriter using the selected headers
			CsvWriter writer = new CsvFileWriter(file, Arrays.asList(fields), true);
			writer.setSeparator(',');
			try
			{
				// write batch to csv
				for (Entity e : entities)
				{
					writer.writeRow(e);
				}
			}
			finally
			{
				writer.close();
			}
		}
	}
}
