import org.apache.commons.io.FileUtils
import org.apache.commons.io.filefilter.WildcardFileFilter
import org.apache.spark.{RangePartitioner, SparkContext, SparkConf, HashPartitioner}
import org.apache.spark.SparkContext._
import sys.process._

import java.io._
import java.nio.file.{Paths, Files}
import java.net._
import java.util.Calendar

import scala.sys.process.Process
import scala.io.Source
import scala.collection.JavaConversions._
import scala.util.Sorting._

import net.sf.samtools._

import tudelft.utils.ChromosomeRange
import tudelft.utils.DictParser
import tudelft.utils.Configuration
import tudelft.utils.SAMRecordIterator

import org.apache.commons.lang3.exception.ExceptionUtils
import org.apache.spark.storage.StorageLevel._

import scala.collection.mutable.{ArrayBuffer, HashMap}

object DNASeqAnalyzer 
{
final val MemString = "-Xmx14336m" 
final val RefFileName = "ucsc.hg19.fasta"
final val SnpFileName = "dbsnp_138.hg19.vcf"
final val ExomeFileName = "gcat_set_025.bed"
//////////////////////////////////////////////////////////////////////////////
def bwaRun (x: String, config: Configuration) :
	Array[(Int, SAMRecord)] =
{
	val refFolder = config.getRefFolder
	val numThreads = config.getNumThreads
	val outFileName = config.getTmpFolder + java.util.UUID.randomUUID().toString.replaceAll("-", "") + ".sam"

	// Create the command string (bwa mem...)and then execute it using the Scala's process package. More help about
	//	Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package.

	// bwa mem refFolder/RefFileName -p -t numOfThreads fastqChunk > outFileName
	val cmdInput = Seq(config.getToolsFolder + "bwa", "mem", refFolder + RefFileName, "-p", "-t", numThreads, x)
	// execute the command
	val exitCode = cmdInput #> new File(outFileName) !

	val bwaKeyValues = new BWAKeyValues(outFileName)
	bwaKeyValues.parseSam()
	val kvPairs: Array[(Int, SAMRecord)] = bwaKeyValues.getKeyValuePairs()

	// Delete the temporary files
	Files.deleteIfExists(Paths.get(outFileName))

	return kvPairs
}

def writeToBAM(fileName: String, samRecordsSorted: Array[SAMRecord], config: Configuration) : ChromosomeRange =
{
	val header = new SAMFileHeader()
	header.setSequenceDictionary(config.getDict())
	val outHeader = header.clone()
	outHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
	val factory = new SAMFileWriterFactory();
	val writer = factory.makeBAMWriter(outHeader, true, new File(fileName));

	val r = new ChromosomeRange()
	val input = new SAMRecordIterator(samRecordsSorted, header, r)
	while(input.hasNext())
	{
		val sam = input.next()
		writer.addAlignment(sam);
	}
	writer.close();

	return r
}

def variantCall (chrRegion: Int, samRecordsSorted: Array[SAMRecord], config: Configuration) :
	Array[(Integer, (Integer, String))] =
{
	val tmpFolder = config.getTmpFolder
	val toolsFolder = config.getToolsFolder
	val refFolder = config.getRefFolder
	val numOfThreads = config.getNumThreads

	// Following is shown how each tool is called. Replace the X in regionX with the chromosome region number (chrRegion).
	// 	You would have to create the command strings (for running jar files) and then execute them using the Scala's process package. More
	// 	help about Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package.
	//	Note that MemString here is -Xmx14336m, and already defined as a constant variable above, and so are reference files' names.

	// SAM records should be sorted by this point
	val regionP1bam = tmpFolder + "region" + chrRegion + "-p1.bam"
	val chrRange = writeToBAM(regionP1bam, samRecordsSorted, config)

	// Picard preprocessing
	//	java MemString -jar toolsFolder/CleanSam.jar INPUT=tmpFolder/regionX-p1.bam OUTPUT=tmpFolder/regionX-p2.bam
	val regionP2bam = tmpFolder + "region" + chrRegion + "-p2.bam"

	Seq("java", MemString, "-jar", toolsFolder + "CleanSam.jar"
		,"INPUT=" + regionP1bam
		,"OUTPUT=" + regionP2bam) !

	//	java MemString -jar toolsFolder/MarkDuplicates.jar INPUT=tmpFolder/regionX-p2.bam OUTPUT=tmpFolder/regionX-p3.bam
	//		METRICS_FILE=tmpFolder/regionX-p3-metrics.txt
	val regionP3bam = tmpFolder + "region" + chrRegion + "-p3.bam"
	val regionP3metrics = tmpFolder + "region" + chrRegion + "-p3-metrics.txt"

	Seq("java", MemString, "-jar", toolsFolder + "MarkDuplicates.jar"
		,"INPUT=" + regionP2bam
		,"OUTPUT=" + regionP3bam
		,"METRICS_FILE=" + regionP3metrics) !

	//	java MemString -jar toolsFolder/AddOrReplaceReadGroups.jar INPUT=tmpFolder/regionX-p3.bam OUTPUT=tmpFolder/regionX.bam
	//		RGID=GROUP1 RGLB=LIB1 RGPL=ILLUMINA RGPU=UNIT1 RGSM=SAMPLE1
	val regionBam = tmpFolder + "region" + chrRegion + ".bam"

	Seq("java", MemString, "-jar", toolsFolder + "AddOrReplaceReadGroups.jar"
		,"INPUT=" + regionP3bam
		,"OUTPUT=" + regionBam
		,"RGID=GROUP1"
		,"RGLB=LIB1"
		,"RGPL=ILLUMINA"
		,"RGPU=UNIT1"
		,"RGSM=SAMPLE1") !

	// 	java MemString -jar toolsFolder/BuildBamIndex.jar INPUT=tmpFolder/regionX.bam

	Seq("java", MemString, "-jar", toolsFolder + "BuildBamIndex.jar", "INPUT=" + regionBam) !

	//	delete tmpFolder/regionX-p1.bam, tmpFolder/regionX-p2.bam, tmpFolder/regionX-p3.bam and tmpFolder/regionX-p3-metrics.txt
	Files.deleteIfExists(Paths.get(regionP1bam))
	Files.deleteIfExists(Paths.get(regionP2bam))
	Files.deleteIfExists(Paths.get(regionP3bam))
	Files.deleteIfExists(Paths.get(regionP3metrics))


	// Make region file
	//	val tmpBed = new File(tmpFolder/tmpX.bed)
	val tmpBed = tmpFolder + "tmp" + chrRegion + ".bed"
	val tmpBedFile = new File(tmpBed)
	chrRange.writeToBedRegionFile(tmpBedFile.getAbsolutePath())

	//	toolsFolder/bedtools intersect -a refFolder/ExomeFileName -b tmpFolder/tmpX.bed -header > tmpFolder/bedX.bed
	val bedBed = tmpFolder + "bed" + chrRegion + ".bed"
	Seq(toolsFolder + "bedtools", "intersect", "-a", refFolder + ExomeFileName, "-b", tmpBed, "-header") #> new File(bedBed) !

	//	delete tmpFolder/tmpX.bed
	Files.deleteIfExists(Paths.get(tmpBed))

	// Indel Realignment
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt numOfThreads -R refFolder/RefFileName
	//		-I tmpFolder/regionX.bam -o tmpFolder/regionX.intervals -L tmpFolder/bedX.bed
	val regionIntervals = tmpFolder + "region" + chrRegion + ".intervals"

	Seq("java", MemString, "-jar", toolsFolder + "GenomeAnalysisTK.jar"
		,"-T", "RealignerTargetCreator", "-nt", numOfThreads, "-R", refFolder + RefFileName
		,"-I", regionBam, "-o", regionIntervals, "-L", bedBed) !

	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T IndelRealigner -R refFolder/RefFileName -I tmpFolder/regionX.bam
	//		-targetIntervals tmpFolder/regionX.intervals -o tmpFolder/regionX-2.bam -L tmpFolder/bedX.bed
	val region2bam = tmpFolder + "region" + chrRegion + "-2.bam"

	Seq("java", MemString, "-jar", toolsFolder + "GenomeAnalysisTK.jar"
			,"-T", "IndelRealigner", "-R", refFolder + RefFileName, "-I", regionBam
			,"-targetIntervals", regionIntervals, "-o", region2bam, "-L", bedBed) !

	//	delete tmpFolder/regionX.bam, tmpFolder/regionX.bai, tmpFolder/regionX.intervals
	Files.deleteIfExists(Paths.get(regionBam))
	Files.deleteIfExists(Paths.get(regionIntervals))
	Files.deleteIfExists(Paths.get(tmpFolder + "region" + chrRegion + ".bai"))

	// Base quality recalibration
	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T BaseRecalibrator -nct numOfThreads -R refFolder/RefFileName -I
	//		tmpFolder/regionX-2.bam -o tmpFolder/regionX.table -L tmpFolder/bedX.bed --disable_auto_index_creation_and_locking_when_reading_rods
	//		-knownSites refFolder/SnpFileName
	val regionTable = tmpFolder + "region" + chrRegion + ".table"
	Seq("java", MemString, "-jar", toolsFolder + "GenomeAnalysisTK.jar"
		,"-T", "BaseRecalibrator" ,"-nct", numOfThreads, "-R", refFolder + RefFileName
		,"-I", region2bam, "-o", regionTable, "-L", bedBed
		,"--disable_auto_index_creation_and_locking_when_reading_rods", "-knownSites", refFolder + SnpFileName) !

	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T PrintReads -R refFolder/RefFileName -I
	//		tmpFolder/regionX-2.bam -o tmpFolder/regionX-3.bam -BSQR tmpFolder/regionX.table -L tmpFolder/bedX.bed
	val region3bam = tmpFolder + "region" + chrRegion + "-3.bam"

	Seq("java", MemString, "-jar", toolsFolder + "GenomeAnalysisTK.jar"
		,"-T", "PrintReads", "-R", refFolder + RefFileName, "-I", region2bam
		,"-o", region3bam, "-BQSR", regionTable, "-L", bedBed) !

	// delete tmpFolder/regionX-2.bam, tmpFolder/regionX-2.bai, tmpFolder/regionX.table
	Files.deleteIfExists(Paths.get(region2bam))
	Files.deleteIfExists(Paths.get(regionTable))
	Files.deleteIfExists(Paths.get(tmpFolder + "region" + chrRegion + "-2.bai"))


	// Haplotype -> Uses the region bed file
	// java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T HaplotypeCaller -nct numOfThreads -R refFolder/RefFileName -I
	//		tmpFolder/regionX-3.bam -o tmpFolder/regionX.vcf  -stand_call_conf 30.0 -stand_emit_conf 30.0 -L tmpFolder/bedX.bed
	//		--no_cmdline_in_header --disable_auto_index_creation_and_locking_when_reading_rods
	val regionVcf = tmpFolder + "region" + chrRegion + ".vcf"

	Seq("java", MemString, "-jar", toolsFolder + "GenomeAnalysisTK.jar"
		,"-T", "HaplotypeCaller", "-nct", numOfThreads, "-R", refFolder + RefFileName
		,"-I", region3bam, "-o", regionVcf, "-stand_call_conf", "30.0", "-stand_emit_conf", "30.0", "-L", bedBed
		,"--no_cmdline_in_header", "--disable_auto_index_creation_and_locking_when_reading_rods") !

	// delete tmpFolder/regionX-3.bam, tmpFolder/regionX-3.bai, tmpFolder/bedX.bed
	Files.deleteIfExists(Paths.get(region3bam))
	Files.deleteIfExists(Paths.get(bedBed))
	Files.deleteIfExists(Paths.get(tmpFolder + "region" + chrRegion + "-3.bai"))

	// return the content of the vcf file produced by the haplotype caller.
	//	Return those in the form of <Chromosome number, <Chromosome Position, line>>
	val vcfParser = new VcfParser(regionVcf)
	vcfParser.parse()
}

def createRegionIndex(numParts: Int, numLines: Int): Array[Int] = {

	val regionIndex = new Array[Int](numLines)

	val div = numLines / numParts
	val remainder = numLines % numParts
	var previousUpperBound = 0

	for (i <- 0 to numParts - 1) {
		if (i < remainder) {
			val lowerBound = previousUpperBound
			val upperBound = (i+1) * div + (i+1)
			previousUpperBound = upperBound
			// update partitionIndex
			for (x <- lowerBound to upperBound -1) {
				regionIndex(x) = i
			}
		} else {
			val lowerBound = previousUpperBound
			val upperBound = (i+1) * div + remainder
			previousUpperBound = upperBound
			// update partitionIndex
			for (x <- lowerBound to upperBound -1) {
				regionIndex(x) = i
			}
		}
	}

	regionIndex
}

def main(args: Array[String]) 
{
	val config = new Configuration()
	config.initialize()
		 
	val conf = new SparkConf().setAppName("DNASeqAnalyzer")
	// For local mode, include the following two lines
	conf.setMaster("local[" + config.getNumInstances() + "]")
	conf.set("spark.cores.max", config.getNumInstances())
	conf.set("spark.ui.port", "9090")
	
	// For cluster mode, include the following commented line
	//conf.set("spark.shuffle.blockTransferService", "nio") 
	
	val sc = new SparkContext(conf)
	
	// Rest of the code goes here

	// --------------------------------
	// PART 1

//	val fastq1 = sc.textFile("../test.txt")
//	val fastq2 = sc.textFile("../test2.txt")

//	val fastq1 = sc.textFile("/data/spark/fastq/fastq1.fq")
//	val fastq2 = sc.textFile("/data/spark/fastq/fastq2.fq")
//
//	val fastq1WithIdx = fastq1.zipWithIndex().map(_.swap)
//	val fastq2WithIdx = fastq2.zipWithIndex().map(_.swap)
//
//	val fastq1Read = fastq1WithIdx.map(x => {
//		val k = x._1.toInt
//		val chunkId = k / 4
//		(chunkId.toLong, x._2)
//	}).groupByKey()
//
//	val fastq2Read = fastq2WithIdx.map(x => {
//		val k = x._1.toInt
//		val chunkId = k / 4
//		(chunkId.toLong, x._2)
//	}).groupByKey()
//
//	val interleave = fastq1Read.union(fastq2Read).reduceByKey(_ ++ _).sortByKey().mapValues(_.mkString("\n"))
//
//	val interleavePart = interleave.partitionBy(new RangePartitioner(config.getNumInstances.toInt, interleave)).cache
//
//	interleavePart.values.saveAsTextFile(config.getInputFolder)

	//--------------------------------------------------
	// PART 2

	//val chunks = sc.wholeTextFiles(config.getInputFolder + "/part-*",config.getNumInstances.toInt)

	// get input filenames
//	val filenames = Seq("ls", config.getInputFolder) #| Seq("grep", "part") !!
//	val fileArray = filenames.split("\n")
//	val chunks = sc.parallelize(fileArray, config.getNumInstances.toInt)
	val chunkDir = new File(config.getInputFolder)
	val chunkFiles = FileUtils.listFiles(chunkDir, new WildcardFileFilter("part*"), null)
	val chunkPaths = chunkFiles.map(f => f.getAbsolutePath).toSeq.sorted
	val chunks = sc.parallelize(chunkPaths, config.getNumInstances.toInt)

	val configBc = sc.broadcast(config)

//	val samRecords = chunks.mapPartitionsWithIndex({
//		(idx, iter) => {
//			iter.map(x => {
//				val receivedConfig = configBc.value
//				val fileName = receivedConfig.getInputFolder + x
//				val bwaOutput = bwaRun(fileName, receivedConfig)
//				bwaOutput
//			})
//		}
//	}, preservesPartitioning = true)

	val samRecords = chunks.mapPartitions({
		iter => {

			var bwaArray = new Array[(Int, SAMRecord)](0)

			val receivedConfig = configBc.value

			for (elem <- iter) {
				val fileName = elem
				val bwaOutput = bwaRun(fileName, receivedConfig)
				bwaArray = bwaArray ++ bwaOutput
			}

			bwaArray.iterator
		}
	}, preservesPartitioning = true)

	//val samRecordsFlat = samRecords.flatMap(x => x)
	val samRecordsFlatWithIdx = samRecords.sortByKey().zipWithIndex().cache()
	val lines = samRecordsFlatWithIdx.count()
	println("lines: " + lines.toString)

	// create region index that will split the samRecord almost equally per instance
	val regionIndex = createRegionIndex(config.getNumInstances.toInt, lines.toInt)
	val regionIndexBc = sc.broadcast(regionIndex)

//	val regionSamRecords = samRecordsFlatWithIdx.map(x => {
//		val receivedIndex = regionIndexBc.value
//		val idx = x._2.toInt
//		(receivedIndex(idx), x._1._2)
//	}).partitionBy(new HashPartitioner(config.getNumInstances.toInt)).cache()

	val regionSamRecords = samRecordsFlatWithIdx.mapPartitions({iter => {
		val receivedIndex = regionIndexBc.value
		iter.map(x => {
			val idx = x._2.toInt
			(receivedIndex(idx), x._1._2)
		})
	}}, preservesPartitioning = true)
		.partitionBy(new HashPartitioner(config.getNumInstances.toInt)).cache()

//	regionSamRecords.saveAsObjectFile(config.getTmpFolder + "regionSamRecords")

//	val regionHist = regionSamRecords.map(x => (x._1, 1)).reduceByKey(_+_)
//	println(regionHist.collect().mkString("\n"))

	val regionSamRecordsGroup = regionSamRecords.groupByKey()

//	val vcf = regionSamRecordsGroup.mapPartitions({
//		iter => {
//			val configReceived = configBc.value
//			// since we are sure that each partition only has one (region, iterable[SAMRecord])
//			// just call iter.next
//			val elem = iter.next()
//			val samArray = elem._2.toArray
//			val region = elem._1
//
//			val samRecordSorted = samArray.sortWith(_.compare(_) < 0)
//			val vcf = variantCall(region, samRecordSorted, configReceived)
//			vcf.iterator
//		}
//	}, preservesPartitioning = true)

	val vcf = regionSamRecordsGroup.mapPartitions({
		iter => {
			val configReceived = configBc.value
			// since we are sure that each partition only has one (region, iterable[SAMRecord])
			// just call iter.next
			var vcfArray = new Array[(Integer, (Integer, String))](0)
			for (elem <- iter) {
				val samArray = elem._2.toArray
				val region = elem._1

				val samRecordSorted = samArray.sortWith(_.compare(_) < 0)
				val vcfResult = variantCall(region, samRecordSorted, configReceived)
				vcfArray = vcfArray ++ vcfResult
			}

			vcfArray.iterator

		}
	}, preservesPartitioning = true).partitionBy(new HashPartitioner(config.getNumInstances.toInt)).cache()

//	val fileArray = Seq("region0.vcf", "region1.vcf"
//		, "region2.vcf", "region3.vcf", "region4.vcf", "region5.vcf", "region6.vcf", "region7.vcf")
//
//	val chunks = sc.parallelize(fileArray, config.getNumInstances.toInt)
//
//	val vcf = chunks.mapPartitions({
//		iter => {
//			val vcfParser = new VcfParser(config.getTmpFolder + iter.next())
//			vcfParser.parse().iterator
//		}
//	}, preservesPartitioning = true)

	// sort the vcf by chromosome number, where in each number the vcf is sorted by position
	val sortedVcf = vcf.groupByKey().mapValues(x => {
		val positionVcfArray = x.toArray
		val positionVcfArraySorted = positionVcfArray.sortWith(_._1 < _._1)
		positionVcfArraySorted
	})

	val sortedVcfFlat = sortedVcf.sortByKey().flatMapValues(x => x)

	// write files
	val vcfLines = sortedVcfFlat.map(x => x._2._2).collect()
	val header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1"

	val outputFile = new File(config.getOutputFolder + "output.vcf")
	outputFile.getParentFile().mkdirs()
	val pw = new PrintWriter(outputFile)

	// write header
	pw.println(header)

	// write vcf
	for (elem <- vcfLines) {
		pw.println(elem)
	}

	pw.close()

}
//////////////////////////////////////////////////////////////////////////////
} // End of Class definition
