import org.apache.spark.Partitioner

/**
 * Created by bagas on 02/11/15.
 */
class EqualSplitPartitioner(numParts: Int, numLines: Int) extends Partitioner {
  override def numPartitions: Int = numParts

  override def getPartition(key: Any): Int = {

    val k = key.asInstanceOf[Long]
    val range = numLines / numParts
    var partitionIdx = 0

    for (i <- 0 to range-1) {
      val lowerBound = i * range
      val upperBound = (i+1) * range
      if ((lowerBound to upperBound-1).contains(k)) {
        partitionIdx = i
      }
    }

    partitionIdx
  }

  override def equals(other: Any): Boolean = other match {
    case esp: EqualSplitPartitioner => esp.numPartitions == numPartitions
    case _ => false
  }
}
