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

import collection.mutable.HashMap

object DNASeqAnalyzer 
{
final val MemString = "-Xmx14336m" 
final val RefFileName = "ucsc.hg19.fasta"
final val SnpFileName = "dbsnp_138.hg19.vcf"
final val ExomeFileName = "gcat_set_025.bed"
//////////////////////////////////////////////////////////////////////////////
//def bwaRun (x: String, config: Configuration) :
//	Array[(Int, SAMRecord)] =
//{
//	val refFolder = config.getRefFolder
//
//	// Create the command string (bwa mem...)and then execute it using the Scala's process package. More help about
//	//	Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package.
//
//	// bwa mem refFolder/RefFileName -p -t numOfThreads fastqChunk > outFileName
//
//	val bwaKeyValues = new BWAKeyValues(outFileName)
//	bwaKeyValues.parseSam()
//	val kvPairs: Array[(Int, SAMRecord)] = bwaKeyValues.getKeyValuePairs()
//
//	// Delete the temporary files
//
//	return kvPairs
//}
//
//def writeToBAM(fileName: String, samRecordsSorted: Array[SAMRecord], config: Configuration) : ChromosomeRange =
//{
//	val header = new SAMFileHeader()
//	header.setSequenceDictionary(config.getDict())
//	val outHeader = header.clone()
//	outHeader.setSortOrder(SAMFileHeader.SortOrder.coordinate);
//	val factory = new SAMFileWriterFactory();
//	val writer = factory.makeBAMWriter(outHeader, true, new File(fileName));
//
//	val r = new ChromosomeRange()
//	val input = new SAMRecordIterator(samRecordsSorted, header, r)
//	while(input.hasNext())
//	{
//		val sam = input.next()
//		writer.addAlignment(sam);
//	}
//	writer.close();
//
//	return r
//}
//
//def variantCall (chrRegion: Int, samRecordsSorted: Array[SAMRecord], config: Configuration) :
//	Array[(Integer, (Integer, String))] =
//{
//	val tmpFolder = config.getTmpFolder
//	val toolsFolder = config.getToolsFolder
//	val refFolder = config.getRefFolder
//	val numOfThreads = config.getNumThreads
//
//	// Following is shown how each tool is called. Replace the X in regionX with the chromosome region number (chrRegion).
//	// 	You would have to create the command strings (for running jar files) and then execute them using the Scala's process package. More
//	// 	help about Scala's process package can be found at http://www.scala-lang.org/api/current/index.html#scala.sys.process.package.
//	//	Note that MemString here is -Xmx14336m, and already defined as a constant variable above, and so are reference files' names.
//
//	// SAM records should be sorted by this point
//	//val chrRange = writeToBAM(tmpFolder/regionX-p1.bam, samRecordsSorted, config)
//
//	// Picard preprocessing
//	//	java MemString -jar toolsFolder/CleanSam.jar INPUT=tmpFolder/regionX-p1.bam OUTPUT=tmpFolder/regionX-p2.bam
//	//	java MemString -jar toolsFolder/MarkDuplicates.jar INPUT=tmpFolder/regionX-p2.bam OUTPUT=tmpFolder/regionX-p3.bam
//	//		METRICS_FILE=tmpFolder/regionX-p3-metrics.txt
//	//	java MemString -jar toolsFolder/AddOrReplaceReadGroups.jar INPUT=tmpFolder/regionX-p3.bam OUTPUT=tmpFolder/regionX.bam
//	//		RGID=GROUP1 RGLB=LIB1 RGPL=ILLUMINA RGPU=UNIT1 RGSM=SAMPLE1
//	// 	java MemString -jar toolsFolder/BuildBamIndex.jar INPUT=tmpFolder/regionX.bam
//	//	delete tmpFolder/regionX-p1.bam, tmpFolder/regionX-p2.bam, tmpFolder/regionX-p3.bam and tmpFolder/regionX-p3-metrics.txt
//
//	// Make region file
//	//	val tmpBed = new File(tmpFolder/tmpX.bed)
//	//	chrRange.writeToBedRegionFile(tmpBed.getAbsolutePath())
//	//	toolsFolder/bedtools intersect -a refFolder/ExomeFileName -b tmpFolder/tmpX.bed -header > tmpFolder/bedX.bed
//	//	delete tmpFolder/tmpX.bed
//
//	// Indel Realignment
//	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt numOfThreads -R refFolder/RefFileName
//	//		-I tmpFolder/regionX.bam -o tmpFolder/regionX.intervals -L tmpFolder/bedX.bed
//	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T IndelRealigner -R refFolder/RefFileName -I tmpFolder/regionX.bam
//	//		-targetIntervals tmpFolder/regionX.intervals -o tmpFolder/regionX-2.bam -L tmpFolder/bedX.bed
//	//	delete tmpFolder/regionX.bam, tmpFolder/regionX.bai, tmpFolder/regionX.intervals
//
//	// Base quality recalibration
//	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T BaseRecalibrator -nct numOfThreads -R refFolder/RefFileName -I
//	//		tmpFolder/regionX-2.bam -o tmpFolder/regionX.table -L tmpFolder/bedX.bed --disable_auto_index_creation_and_locking_when_reading_rods
//	//		-knownSites refFolder/SnpFileName
//	//	java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T PrintReads -R refFolder/RefFileName -I
//	//		tmpFolder/regionX-2.bam -o tmpFolder/regionX-3.bam -BSQR tmpFolder/regionX.table -L tmpFolder/bedX.bed
//	// delete tmpFolder/regionX-2.bam, tmpFolder/regionX-2.bai, tmpFolder/regionX.table
//
//	// Haplotype -> Uses the region bed file
//	// java MemString -jar toolsFolder/GenomeAnalysisTK.jar -T HaplotypeCaller -nct numOfThreads -R refFolder/RefFileName -I
//	//		tmpFolder/regionX-3.bam -o tmpFolder/regionX.vcf  -stand_call_conf 30.0 -stand_emit_conf 30.0 -L tmpFolder/bedX.bed
//	//		--no_cmdline_in_header --disable_auto_index_creation_and_locking_when_reading_rods
//	// delete tmpFolder/regionX-3.bam, tmpFolder/regionX-3.bai, tmpFolder/bedX.bed
//
//	// return the content of the vcf file produced by the haplotype caller.
//	//	Return those in the form of <Chromsome number, <Chromosome Position, line>>
//}

def main(args: Array[String]) 
{
	val config = new Configuration()
	config.initialize()
		 
	val conf = new SparkConf().setAppName("DNASeqAnalyzer")
	// For local mode, include the following two lines
	conf.setMaster("local[" + config.getNumInstances() + "]")
	conf.set("spark.cores.max", config.getNumInstances())
	
	// For cluster mode, include the following commented line
	//conf.set("spark.shuffle.blockTransferService", "nio") 
	
	val sc = new SparkContext(conf)
	
	// Rest of the code goes here

	val fastq1 = sc.textFile("../test.txt")
	val fastq2 = sc.textFile("../test2.txt")
	val numLines = fastq1.count().toInt

	val fastq1WithIdx = fastq1.zipWithIndex().map(_.swap)
	val fastq2WithIdx = fastq2.zipWithIndex().map(_.swap)

	val fastq1Read = fastq1WithIdx.map(x => {
		val k = x._1.toInt
		val range = numLines / 4
		var chunkId = 0
		for (i <- 0 to range - 1) {
			val lowerBound = i * range
			val upperBound = (i+1) * range
			if ((lowerBound to upperBound -1).contains(k)) chunkId = i
		}
		(chunkId.toLong, x._2)
	}).groupByKey()

	val fastq2Read = fastq2WithIdx.map(x => {
		val k = x._1.toInt
		val range = numLines / 4
		var chunkId = 0
		for (i <- 0 to range - 1) {
			val lowerBound = i * range
			val upperBound = (i+1) * range
			if ((lowerBound to upperBound -1).contains(k)) chunkId = i
		}
		(chunkId.toLong, x._2)
	}).groupByKey()

	val interleave = fastq1Read.union(fastq2Read).reduceByKey(_ ++ _)

	val interleavePart = interleave.partitionBy(new RangePartitioner(config.getNumInstances.toInt, interleave)).cache

	interleavePart.saveAsTextFile("output")

//	val mapper = interleavePart.mapPartitionsWithIndex({
//		(idx, iter) => {
//			val pw = new PrintWriter(idx.toString + ".chunk")
//			for (x <- iter) {
//				val list = x._2
//				for (elem <- list) pw.println(elem)
//			}
//			iter
//		}
//	}, preservesPartitioning = true)
//
//	mapper.take(100).foreach(println)


}
//////////////////////////////////////////////////////////////////////////////
} // End of Class definition
