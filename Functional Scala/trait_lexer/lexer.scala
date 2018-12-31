
class Pos private(val prog: String, val offs: Int, val line: Int, val col: Int) {
  def this(prog: String) = this(prog, 0, 1, 1)

  def ch = if (offs == prog.length) -1 else prog(offs)

  def inc = ch match {
    case '\n' => new Pos(prog, offs+1, line+1, 1)
    case -1   => this
    case _    => new Pos(prog, offs+1, line, col+1)
  }

  override def toString = "(" + line + ", " + col + ")"
}

object DomainTags extends Enumeration {
  type Tag = Value
  val WHITESPACE, IDENT, NUMBER, END_OF_PROGRAM = Value
}
import DomainTags._

class Token(val start: Pos, scanner: Scanner) {
  val (tag, follow) = start.ch match {
    case -1 => (END_OF_PROGRAM, start)
    case _  => scanner.scan(start)
  }

  def image = start.prog.substring(start.offs, follow.offs)

  def next = new Token(follow, scanner)
}


class Scanner {
  def scan(start: Pos): (Tag, Pos) =
    sys.error("syntax error at " + start)
}

trait Whitespaces extends Scanner {
  private def missWhitespace(pos: Pos): Pos = pos.ch match {
    case ' '  => missWhitespace(pos.inc)
    case '\t' => missWhitespace(pos.inc)
    case '\n' => missWhitespace(pos.inc)
    case _    => pos
  }

  override def scan(start: Pos) = {
    val follow = missWhitespace(start)
    if (start != follow) (WHITESPACE, follow)
    else super.scan(start)
  }
}

var t = new Token(
  new Pos("1234\n   alpha"),
  new Scanner with Numbers with Whitespaces
)

while (t.tag != END_OF_PROGRAM) {
  println(t.tag.toString + " " + t.start + "-" + t.follow + ": " + t.image)
  t = t.next
}
