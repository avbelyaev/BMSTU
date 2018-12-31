#include <stdint.h>
#include "llvm/ADT/STLExtras.h"
#include "llvm/IR/IRBuilder.h"
#include "llvm/IR/Module.h"
#include "llvm/IR/Verifier.h"

using namespace std;

#define TOK_EOF 1
#define TOK_IDENT 2
#define TOK_NUM 3
#define TOK_IF 6
#define TOK_ELSE 7

static string IdentifierStr;
static int NumVal;
static int CurToken;

static int GetToken() {

    static int posCh = ' ';

    while (isspace(posCh))
        posCh = getchar();
    if (isalpha(posCh)) {
        IdentifierStr = posCh;
        while (isalnum((posCh = getchar())))
            IdentifierStr += posCh;

        if (IdentifierStr == "if")
            return TOK_IF;
        if (IdentifierStr == "else")
            return TOK_ELSE;

        return TOK_IDENT;
    }

    if (isdigit(posCh)) {
        string NumStr;
        do {
            NumStr += posCh;
            posCh = getchar();
        } while (isdigit(posCh));

        NumVal = strtod(NumStr.c_str(), nullptr);
        return TOK_NUM;
    }

    if (posCh == EOF)
        return TOK_EOF;

    int thisChar = posCh;
    posCh = getchar();

    return thisChar;
}

void NextToken(){
    CurToken = GetToken();
}




//===-----------------------------------------------------===//
//=---                       EXPR                        ---=//
//===-----------------------------------------------------===//
class ExprAST {
public:
    virtual ~ExprAST() {}
    virtual llvm::Value *codegen() = 0;
};



//===-----------------------------------------------------===//
//=---                     NUMBER                        ---=//
//===-----------------------------------------------------===//
class NumberExprAST : public ExprAST {
    int Val;

public:
    NumberExprAST(int Val) : Val(Val) {}
    llvm::Value *codegen() override;
};


//===-----------------------------------------------------===//
//=---                     ASSIGN                        ---=//
//===-----------------------------------------------------===//
class AssignAST : public ExprAST {
    string Name;
    ExprAST *assExpr;
public:
    AssignAST(const string &Name, ExprAST *Expr) {
        this->Name = Name;
        this->assExpr = Expr;
    }
    llvm::Value *codegen() override;
};


//===-----------------------------------------------------===//
//=---                       VAR                         ---=//
//===-----------------------------------------------------===//
class VariableExprAST : public ExprAST {
    string Name;
public:
    VariableExprAST(const string &Name) : Name(Name) {}
    llvm::Value *codegen() override;
};


//===-----------------------------------------------------===//
//=---                        IF                         ---=//
//===-----------------------------------------------------===//
class IfExprAST : public ExprAST {
    ExprAST *Cond, *Then, *Else;
public:
    IfExprAST(ExprAST *Cond, ExprAST *Then, ExprAST *Else){
        this->Cond = Cond;
        this->Then = Then;
        this->Else = Else;
    }
    llvm::Value *codegen() override;
};


//===-----------------------------------------------------===//
//=---                     BINARY                        ---=//
//===-----------------------------------------------------===//
class BinaryExprAST : public ExprAST {
    char Op;
    ExprAST *LHS, *RHS;

public:
    BinaryExprAST(char Op, ExprAST *LHS,
                  ExprAST *RHS)
            : Op(Op), LHS(LHS), RHS(RHS) {}
    llvm::Value *codegen() override;
};





//===---------------------========------------------------===//
//=---                      PARSE                        ---=//
//===---------------------========------------------------===//

//создание Абстр Синт Дерева

//https://habrahabr.ru/post/120005/

static ExprAST *ParseExpression();

static ExprAST *ParseNumberExpr() {
    NumberExprAST *Result = new NumberExprAST(NumVal);
    NextToken();
    return Result;
}

static ExprAST *ParseParenExpr() {
    NextToken();
    ExprAST *V = ParseExpression();
    if (!V)
        return nullptr;

    if (CurToken != ')') {
        printf("%s\n", "expected ')'");
        return nullptr;
    }
    NextToken();
    return V;
}

static ExprAST *ParseVarExpr() {
    string Name = IdentifierStr;

    NextToken();

    if (CurToken != '=') {
        return new VariableExprAST(Name);
    }

    NextToken();

    ExprAST *varDecl = ParseExpression();
    if (!varDecl) {
        printf("expected expr after \"=\"!");
        return nullptr;
    }
    return new AssignAST(Name, varDecl);
}

static ExprAST *ParseIfExpr() {
    NextToken();
    if (CurToken != '(') {
        printf("%s %c %d\n,", "expected '(' after if, no ", CurToken, CurToken);
    }
    NextToken();
    ExprAST *Cond = ParseExpression();
    if (!Cond)
        return nullptr;

    if (CurToken != ')') {
        printf("%s %c %d\n,", "expected ')' after if", CurToken, CurToken);
    }

    NextToken();
    NextToken();

    ExprAST *Then = ParseExpression();
    NextToken();
    ExprAST *Else = ParseExpression();

    return new IfExprAST(Cond, Then, Else);
}

static ExprAST *ParsePrimary() {
    switch (CurToken) {
        case TOK_IDENT:
            return ParseVarExpr();
        case TOK_NUM:
            return ParseNumberExpr();
        case '(':
            return ParseParenExpr();
        case TOK_IF:
            return ParseIfExpr();
        default:
            printf("unknown token when expecting an expression: %c\n", CurToken);
            return nullptr;
    }
}


static ExprAST *ParseBinary (ExprAST *expr) {
    ExprAST *LHS = expr;
    while (1) {
        //  printf("parsemath1c: %d",CurToken);
        if (CurToken!='+' && CurToken!='-')
            return LHS;
        char Op = CurToken;
        NextToken();
        //printf("parsemath2c: %d",CurToken);
        ExprAST *RHS = ParsePrimary();
        if (!RHS)
            return nullptr;
        LHS = new BinaryExprAST(Op, LHS, RHS);
    }
}

static ExprAST *ParseExpression() {
    ExprAST *expr = ParsePrimary();
    if (!expr) return nullptr;
    return ParseBinary(expr);
}





//===---------------------========------------------------===//
//=---                    GENERATE                       ---=//
//===---------------------========------------------------===//

//https://habrahabr.ru/post/120424/

//генерация Intermediate Representation по AST:

//определим виртуальные методы кодогенерации (Codegen) для каждого AST-класса
// Метод Codegen() вернёт IR для данного узла AST вместе со всеми зависимыми от него, и все они возвращают объект LLVM Value.
// "Value" является классом, используемым для представления «регистра Static Single Assigment (SSA)» или «Значения SSA» в LLVM.

//TheModule является конструкцией, содержащей все функции и глобальные переменные в куске кода.
// В большинстве случаев, это структуры верхнего уровня, которые использует LLVM IR для содержащегося кода.

//Объект Builder является вспомогательным объектом, который позволяет генерировать инструкции LLVM

//Карта NamedValues отслеживает, какие значения определены в текущей области видимости и каково их LLVM-представление.
// (Другими словами, это таблица символов для кода). В текущей форме в Kaleidoscope, единственное, что может ссылаться —
// это параметры функций. Таким образом, в этой карте при генерации кода для тела функции будут расположены параметры этой функции.



static llvm::Module *MainModule;
static llvm::Function *MainFunction;
static llvm::LLVMContext MainContext;
static llvm::IRBuilder<> Builder(llvm::getGlobalContext());
static map<string, llvm::Value *> NamedValues;

llvm::Value *AssignAST::codegen(){
    llvm::Value *v = assExpr->codegen();
    if (!v)
        return nullptr;

    llvm::AllocaInst *Alloca = Builder.CreateAlloca(llvm::Type::getInt32Ty(llvm::getGlobalContext()), 0, Name.c_str());
    Builder.CreateStore(v, Alloca);
    llvm::Value *CurVar = Builder.CreateLoad(Alloca);

    NamedValues[Name] = Alloca;
    return CurVar;
}

llvm::Value *NumberExprAST::codegen() {
    return llvm::ConstantInt::get(llvm::getGlobalContext(), llvm::APInt(32, Val, false));
    //APInt(unsigned numBits, uint64_t val, bool isSigned = false)
}

llvm::Value *VariableExprAST::codegen() {
    llvm::Value *V = NamedValues[Name];
    if (!V) {
        printf("unknown variable name\n");
        return nullptr;
    }
    return Builder.CreateLoad(V, Name);
}

llvm::Value *BinaryExprAST::codegen() {
    llvm::Value *L = LHS->codegen();
    llvm::Value *R = RHS->codegen();
    if (!L || !R)
        return nullptr;
    switch (Op) {
        case '+':
            return Builder.CreateAdd(L, R, "addtmp");
        case '-':
            return Builder.CreateSub(L, R, "subtmp");
        default:
            printf("%s\n","invalid binary operator");
            return nullptr;
    }
    //IRBuilder знает, куда вставлять вновь созданные инструкции, и всё, что нужно сделать вам —
    // это указать, какие инструкции создавать (например, "CreateFAdd"),
    // какие операнды использовать (в данном случае L и R) и, если требуется, то какое использовать имя для генерируемой инструкции.
}

llvm::Value *IfExprAST::codegen() {
    llvm::Value *CondV = Cond->codegen();
    if (!CondV)
        return nullptr;

    //https://habrahabr.ru/post/120881/
    // Конвертация условия в булево сравнением с 0.0.
    CondV = Builder.CreateICmpNE(
            CondV, llvm::ConstantInt::get(llvm::getGlobalContext(), llvm::APInt(32, 0, false)), "ifcond");

    //Первая строка получает текущий объект Function (формируемой функции).
    // Она получает его, спрашивая Builder о текущем базовом блоке и получая его «родителя» (текущую функцию).

    //Затем он создает три блока. Обратите внимание, что он передаёт "TheFunction" в конструктор блока "then".
    // Это вызывает конструктор с автоматической вставкой нового блока в конец указанной функции.
    // Два других блока тоже создаются, но еще не вставлены в функцию.
    llvm::Function *MainFunction = Builder.GetInsertBlock()->getParent();

    llvm::BasicBlock *ThenBB = llvm::BasicBlock::Create(llvm::getGlobalContext(), "then", MainFunction);
    llvm::BasicBlock *ElseBB = llvm::BasicBlock::Create(llvm::getGlobalContext(), "else");
    llvm::BasicBlock *MergeBB = llvm::BasicBlock::Create(llvm::getGlobalContext(), "ifcont");

    Builder.CreateCondBr(CondV, ThenBB, ElseBB);

    // Генерируем значение.
    Builder.SetInsertPoint(ThenBB);
    llvm::Value *ThenV = Then->codegen();
    if (!ThenV) return nullptr;

    Builder.CreateBr(MergeBB);
    // Кодогенерация 'Then' может изменить текущий блок, обновляем ThenBB для PHI.
    ThenBB = Builder.GetInsertBlock();

    // Генерируем блок else.
    MainFunction->getBasicBlockList().push_back(ElseBB);
    Builder.SetInsertPoint(ElseBB);
    llvm::Value *ElseV = Else->codegen();
    if (!ElseV) return nullptr;

    Builder.CreateBr(MergeBB);
    // Кодогенерация 'Else' может изменить текущий блок, обновляем ElseBB для PHI.
    ElseBB = Builder.GetInsertBlock();


    // Генерация блока слияния.
    MainFunction->getBasicBlockList().push_back(MergeBB);
    Builder.SetInsertPoint(MergeBB);
    llvm::PHINode *PN = Builder.CreatePHI(llvm::Type::getInt32Ty(llvm::getGlobalContext()), 2, "iftmp");

    PN->addIncoming(ThenV, ThenBB);
    PN->addIncoming(ElseV, ElseBB);
    return PN;
}

static llvm::Value *Parse() {
    NextToken();

    ExprAST *expr = nullptr;
    llvm::Value *RetVal = nullptr;
    while (CurToken != TOK_EOF) {
        expr = ParseExpression();
        if (expr) {
            RetVal = expr->codegen();
            if (!RetVal) NextToken();
        } else {
            return nullptr;
        }
    }

    return RetVal;
}


int main() {
    // Создаём модуль, который будет хранить весь код.
    MainModule = new llvm::Module("main", MainContext);

    MainFunction = MainModule->getFunction("main");

    llvm::FunctionType *FT =
            llvm::FunctionType::get(llvm::Type::getVoidTy(llvm::getGlobalContext()),false);
    MainFunction =
            llvm::Function::Create(FT, llvm::GlobalValue::CommonLinkage, "main", MainModule);

    llvm::BasicBlock *BB = llvm::BasicBlock::Create(MainContext, "entry", MainFunction);
    Builder.SetInsertPoint(BB);

    llvm::Value *ret = Parse();
    Builder.CreateRet(ret);
    //llvm::verifyFunction(*MainFunction);

    // Выводим сгенерированный код.
    MainModule->dump();

    return 0;
}


//Локальные значения обозначаются префиксом %, а глобальные — @.
// Локальные значения также называют регистрами, а LLVM — виртуальной машиной с бесконечным числом регистров.

//что из себя представляет LLVM IR.
// В одном предложении его можно охарактеризовать как типизированный трёхадресный код в SSA-форме.

//базовые блоки, содержат последовательность инструкций, заканчивающуюся инструкцией-терминатором,
// ЯВНО передающей управление в другой блок.
// Базовые блоки в LLVM обозначаются метками, а терминаторами являются следующие инструкции: br, ret, switch ...

//Обратиться к памяти можно только с помощью двух инструкций, названия которых говорят сами за себя: load и store
// Выделение памяти - alloca. Память выделяется на стеке. (malloc - в куче)
// Память, выделенная alloca, автоматически освобождается при выходе из функции при помощи инструкций ret или unwind.
// http://llvm.org/docs/LangRef.html


//g++ -g -O3 main.cpp `llvm-config --cxxflags --ldflags --system-libs --libs core` -fno-rtti -o prog
//./prog

/*
 * a = 2
b = 3
if (a) then
b = 4
else
b = a + 5
d = b
; ModuleID = 'main'

define common void @main() {
entry:
  %a = alloca i32
  store i32 2, i32* %a
  %0 = load i32, i32* %a
  %b = alloca i32
  store i32 3, i32* %b
  %1 = load i32, i32* %b
  %a1 = load i32, i32* %a
  %ifcond = icmp ne i32 %a1, 0
  br i1 %ifcond, label %then, label %else

then:                                             ; preds = %entry
  %b2 = alloca i32
  store i32 4, i32* %b2
  %2 = load i32, i32* %b2
  br label %ifcont

else:                                             ; preds = %entry
  %a3 = load i32, i32* %a
  %addtmp = add i32 %a3, 5
  %b4 = alloca i32
  store i32 %addtmp, i32* %b4
  %3 = load i32, i32* %b4
  br label %ifcont

ifcont:                                           ; preds = %else, %then
  %iftmp = phi i32 [ %2, %then ], [ %3, %else ]
  %b5 = load i32, i32* %b4
  %d = alloca i32
  store i32 %b5, i32* %d
  %4 = load i32, i32* %d
  ret i32 %4
}
 */

//CTRL+D