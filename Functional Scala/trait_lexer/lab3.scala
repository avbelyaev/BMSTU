
class Pos private(val prog: String, val offs: Int, val line: Int, val col: Int) {
	def this(prog: String) = this(prog, 0, 1, 1)

	def ch = if (offs == prog.length) -1 else prog.codePointAt(offs)

	def inc = ch match {
		case '\n' => new Pos(prog, offs+1, line+1, 1)
		case -1   => this
		case _    => new Pos(prog, offs+1, line, col+1)
	}

	override def toString = "(" + line + ", " + col + ")"

	def isLetter = if (ch >= 'a' && ch <= 'z') true else false

	def isDigit = if (ch >= '0' && ch <= '9') true else false
}

object DomainTags extends Enumeration {
	type Tag = Value
	val WHITESPACE, IDENT, NUMBER, OPERATION, ERROR, END_OF_PROGRAM = Value
}

import DomainTags._

class Scanner {
	def scan(start: Pos): (Tag, Pos) = {
		//sys.error("Syntax error at " + start)
		 (ERROR, start.inc)
	}
}

class Token(val start: Pos, scanner: Scanner) {
	val (tag, follow) = start.ch match {
	  case -1 => (END_OF_PROGRAM, start)
	  case _  => scanner.scan(start)
	}

	def image = start.prog.substring(start.offs, follow.offs)

	def next = new Token(follow, scanner)
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


// Вар 7
// Идентификаторы: последовательности латинских букв и десятичных цифр, оканчивающиеся на цифру.
// Числовые литералы: последовательности десятичных цифр, органиченные знаками «<» и «>».
// Операции: «<=», «=», «==».


trait Idents extends Scanner {

	private def recognizeIdent(pos: Pos): (Pos, Boolean) = pos.ch match {
		case _ if pos.isDigit &&
		 (!pos.inc.isLetter && !pos.inc.isDigit) 	=> (pos.inc, false)
		case _ if pos.isDigit                     => recognizeIdent(pos.inc)
		case _ if pos.isLetter							=> recognizeIdent(pos.inc)
		case _ 												=> (pos, true)
	}

	override def scan(start: Pos) = {
		val (follow, err) =
			if (start.isDigit || start.isLetter) recognizeIdent(start)
			else (start, false)

		if (start != follow) {
			if (err) (ERROR, follow)
			else (IDENT, follow)
		}
		else super.scan(start)
	}
}



trait Numbers extends Scanner {

	private def recognizeNumber(pos: Pos): (Pos, Boolean) = pos.ch match {
		case '<'     				=> recognizeNumber(pos.inc)
		case '>'     				=> (pos.inc, false)
		case _ if pos.isDigit 	=> recognizeNumber(pos.inc)
		case _ 						=> (pos, true)
	}

	override def scan(start: Pos) = {
		val (follow, err) =
			if ('<' == start.ch && start.inc.isDigit) recognizeNumber(start)
			else (start, false)

		if (start != follow) {
			if (err) (ERROR, follow)
			else (NUMBER, follow)
		}
		else super.scan(start)
	}
}



trait Operations extends Scanner {

	override def scan(start: Pos) = {
		val follow =
			if ('<' == start.ch && '=' == start.inc.ch) start.inc.inc
			else if ('=' == start.ch&& '=' == start.inc.ch) start.inc.inc
				else if ('=' == start.ch) start.inc
					else start

		if (start != follow) (OPERATION,  follow)
		else super.scan(start)
  }
}

//X extends A with B, C, D.
//L(X) = X, D, C, B, A
var t = new Token(
	new Pos(" <   <456 <== <1337> == 4w3s0m3 1den7> 3rr0r "),
	new Scanner
	  with Numbers
	  with Operations
	  with Idents
	  with Whitespaces
)

while (t.tag != END_OF_PROGRAM) {
	println(t.tag.toString + " " + t.start + "-" + t.follow + ": \t\t" + t.image)
	t = t.next
}