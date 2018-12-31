package main

import (
	"fmt"
	"go/ast"
	"go/format"
	"go/parser"
	"go/token"
	"io/ioutil"
	"os"
	"sort"
)

type storage struct {
	_type *ast.FuncType
	_body *ast.BlockStmt
}

//=====================================================
//-----------------GLOBAL-FUNCTION---------------------
//=====================================================

func insertFunc(file *ast.File, name string, intrails storage) {

	//compose function of intrails and name
	new_func := &ast.FuncDecl{
		Doc:  nil,
		Recv: nil,
		Name: &ast.Ident{
			NamePos: token.NoPos,
			Name:    name,
		},
		Type: intrails._type,
		Body: intrails._body,
	}

	//add it to the back
	file.Decls = append(file.Decls, new_func)
}

//=====================================================
//----------------------IMPORT-------------------------
//=====================================================

func insertImport(file *ast.File, import_name string) {

	//import was previusly deleted by ast.FilterFile()
	new_impt := &ast.GenDecl{
		Doc:    nil,
		TokPos: token.NoPos,
		Tok:    token.IMPORT,
		Lparen: token.NoPos,
		Specs: []ast.Spec{
			&ast.ImportSpec{
				Doc:  nil,
				Name: nil,
				Path: &ast.BasicLit{
					Kind:  token.STRING,
					Value: "\"" + import_name + "\"",
				},
				Comment: nil,
			},
		},
		Rparen: token.NoPos,
	}

	var before, after []ast.Decl

	if len(file.Decls) > 0 {
		hasImport := false
		if genDecl, ok := file.Decls[0].(*ast.GenDecl); ok {
			hasImport = genDecl.Tok == token.IMPORT
		}

		if hasImport {
			before, after = []ast.Decl{file.Decls[0]}, file.Decls[1:]
		} else {
			after = file.Decls
		}
	}

	//post it to the beginning
	file.Decls = append(before, new_impt)
	file.Decls = append(file.Decls, after...)
}

//=====================================================
//---------------------INT-var-------------------------
//=====================================================

func insertIntVar(file *ast.File, name string, value int) {
	var before, after []ast.Decl

	if len(file.Decls) > 0 {
		hasImport := false
		if genDecl, ok := file.Decls[0].(*ast.GenDecl); ok {
			hasImport = genDecl.Tok == token.IMPORT
		}

		if hasImport {
			before, after = []ast.Decl{file.Decls[0]}, file.Decls[1:]
		} else {
			after = file.Decls
		}
	}

	file.Decls = append(before,
		&ast.GenDecl{
			Tok: token.VAR,
			Specs: []ast.Spec{
				&ast.ValueSpec{
					Names: []*ast.Ident{ast.NewIdent(name)},
					Type:  ast.NewIdent("int"),
					Values: []ast.Expr{
						&ast.BasicLit{
							Kind:  token.INT,
							Value: fmt.Sprintf("%d", value),
						},
					},
				},
			},
		},
	)

	file.Decls = append(file.Decls, after...)

}

//=====================================================
//-----------------------CONST-------------------------
//=====================================================

func insertConst(file *ast.File, name string, value string) {

	var before, after []ast.Decl

	if len(file.Decls) > 0 {
		hasImport := false
		if genDecl, ok := file.Decls[0].(*ast.GenDecl); ok {
			hasImport = genDecl.Tok == token.IMPORT
		}

		if hasImport {
			before, after = []ast.Decl{file.Decls[0]}, file.Decls[1:]
		} else {
			after = file.Decls
		}
	}

	file.Decls = append(before,
		&ast.GenDecl{
			Tok: token.CONST,
			Specs: []ast.Spec{
				&ast.ValueSpec{
					Doc:   nil,
					Names: []*ast.Ident{ast.NewIdent(name)},
					Type:  nil,
					Values: []ast.Expr{
						&ast.BasicLit{
							Kind:  token.STRING,
							Value: "\"" + value + "\"",
						},
					},
					Comment: nil,
				},
			},
		},
	)

	file.Decls = append(file.Decls, after...)
}

//=====================================================
//------------------Printf(hello)----------------------
//=====================================================

func insertHello(file *ast.File) {
	ast.Inspect(file, func(node ast.Node) bool {
		if ifStmt, ok := node.(*ast.IfStmt); ok {
			ifStmt.Body.List = append(
				[]ast.Stmt{
					&ast.ExprStmt{
						X: &ast.CallExpr{
							Fun: &ast.SelectorExpr{
								X:   ast.NewIdent("fmt"),
								Sel: ast.NewIdent("Printf"),
							},
							Args: []ast.Expr{
								&ast.BasicLit{
									Kind:  token.CONST,
									Value: "\"hello\"",
								},
							},
						},
					},
				},
				ifStmt.Body.List...,
			)
		}
		return true
	})
}

func modify_tree(file *ast.File, fs *token.FileSet, src []byte) {

	str := string(src)
	m_cont := make(map[string]storage) //functions map [func_name]<func_body as container>

	ast.Inspect(file, func(node ast.Node) bool {

		switch x := node.(type) {

		case *ast.FuncDecl:

			//Вар 6
			//Объявления глобальных функций должны располагаться
			//в конце программы в алфавитном порядке

			func_name := x.Name.Name
			func_itself := x.Body
			func_type := x.Type

			func_start_offs := fs.Position(func_itself.Pos()).Offset
			func_end_offs := fs.Position(func_itself.End()).Offset
			func_line := fs.Position(func_itself.Pos()).Line
			func_filename := fs.Position(func_itself.Pos()).Filename
			func_body := str[func_start_offs:func_end_offs]

			if "main" != func_name {

				t := storage{func_type, func_itself}
				m_cont[func_name] = t

				fmt.Printf("FuncName:[%s]\nFileName:[%s]\nLine:[%d] Offset:[%d..%d]\n FuncBody:[%s]\n",
					func_name, func_filename, func_line, func_start_offs, func_end_offs, func_body)
			}

		}

		return true
	})

	//remove all nodes that were previously saved to map
	for key, _ := range m_cont {

		remove_node(file, key)
	}

	//sort func_names alphabetically
	sorted := make([]string, len(m_cont))
	i := 0

	for k, _ := range m_cont {
		sorted[i] = k
		i++
	}
	sort.Strings(sorted)

	//add functions to the back
	for i := 0; i < len(m_cont); i++ {

		name := sorted[i]
		body := m_cont[name]
		insertFunc(file, name, body)

		fmt.Printf("Func [%s] were added\n", name)
	}

	//restore Import statemnt that was deleted by ast.FilterFile()
	insertImport(file, "fmt")

}

//=====================================================
//--------------------REMOVE-NODE----------------------
//=====================================================

func remove_node(file *ast.File, name string) {
	//as said in filter.go:
	// FilterFile trims the AST for a Go file in place by removing all
	// names from top-level declarations that don't pass through the filter f...
	// Import declarations are always removed.
	ast.FilterFile(file, func(ident string) bool {
		println("Testing filter for: ", ident)
		//keep all elements that dont match the name
		return ident != name
	})
}

func main() {
	if len(os.Args) != 2 {
		fmt.Printf("usage: demo.exe <filename.go>\n")
		return
	}

	// Создаём хранилище данных об исходных файлах
	fset := token.NewFileSet()
	// Сохраняем исходный файл как строку
	src, err := ioutil.ReadFile(os.Args[1])

	if nil != err {
		fmt.Printf("error reading file!\n")
	}

	// Вызываем парсер
	if file, err := parser.ParseFile(
		fset,                 // данные об исходниках
		os.Args[1],           // имя файла с исходником программы
		src,                  // пусть парсер сам загрузит исходник
		parser.ParseComments, // приказываем сохранять комментарии
	); err == nil {
		// Если парсер отработал без ошибок, печатаем дерево
		//insertIntVar(file, "xxx", 666)
		//insertHello(file)
		modify_tree(file, fset, src)

		if format.Node(os.Stdout, fset, file) != nil {
			fmt.Printf("Formatter error: %v\n", err)
		}

		//ast.Fprint(os.Stdout, fset, file, nil)

		// Собираем файл из модифицированного дерева
		tmp := string(os.Args[1])
		f, _ := os.Create(tmp[:len(tmp)-3] + "_MOD.go")
		format.Node(f, fset, file)

		fmt.Printf("%s_MOD.go has been created\n", tmp[:len(tmp)-3])

	} else {
		// в противном случае, выводим сообщение об ошибке
		fmt.Printf("Error: %v", err)
	}
}
