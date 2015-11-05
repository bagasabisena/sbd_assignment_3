import org.apache.spark.Partitioner

/**
 * Created by bagas on 02/11/15.
 */
class EqualSplitPartitioner(numParts: Int, numLines: Int) extends Partitioner {

  val partitionIndex = createPartitionIndex(numParts, numLines)

  def createPartitionIndex(numParts: Int, numLines: Int): Array[Int] = {

    val partitionIndex = new Array[Int](numLines)

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
          partitionIndex(x) = i
        }
      } else {
        val lowerBound = previousUpperBound
        val upperBound = (i+1) * div + remainder
        previousUpperBound = upperBound
        // update partitionIndex
        for (x <- lowerBound to upperBound -1) {
          partitionIndex(x) = i
        }
      }
    }

    partitionIndex
  }

  override def numPartitions: Int = numParts

  override def getPartition(key: Any): Int = {

    val k = key.asInstanceOf[Long]
//    val divider = numLines / numParts
//    val modulo = numLines % numParts
//
//    val partitionIdx = (k-modulo) / divider
//    partitionIdx.toInt
    partitionIndex(k.toInt)
  }

  override def equals(other: Any): Boolean = other match {
    case esp: EqualSplitPartitioner => esp.numPartitions == numPartitions
    case _ => false
  }
}
