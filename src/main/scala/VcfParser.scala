import scala.collection.mutable.ArrayBuffer
import scala.io.Source

/**
 * Created by bagas on 06/11/15.
 */
class VcfParser(filename:String) {

  val vcfIterator = Source.fromFile(filename).getLines()
  var header:String = _

  def parse(): Array[(Integer, (Integer, String))] = {

    // ignore ## lines
    var elem: String = vcfIterator.next()
    while (elem.substring(0,2) == "##") {
      elem = vcfIterator.next()
    }

    // get header
    header = elem

    // parse the rest
    val arr = new ArrayBuffer[(Integer, (Integer, String))]()
    for (elem <- vcfIterator) {
      val line = elem
      val elemArray = line.split("\t")
      val chrNumber = getChromosomeNumber(elemArray(0))
      val position = elemArray(1).toInt
      arr += new Pair(chrNumber, new Pair(position, line))
    }

    arr.toArray
  }

  def getChromosomeNumber(chrString: String): Int = {

    // Strip "chr" from chrX
    val numString = chrString.substring(3)
    val chrNumber:Int = numString match {
      case "X" => 23
      case "Y" => 24
      case _ => numString.toInt
    }

    chrNumber

  }

}
