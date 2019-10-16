-- noinspection SqlNoDataSourceInspectionForFile

-- ### 3. Создать tablespace
-- 1. Создать tablespace “myts1” (myts2) в директории MyDB1 (MyDB2) с пользователем postgres
/*
- в контейнере:
mkdir -p /home/space1 && chown postgres /home/space1
mkdir -p /home/space2 && chown postgres /home/space2
*/
-- - по отдельности запускаем:
CREATE TABLESPACE myts1 OWNER ics LOCATION '/home/space1';
CREATE TABLESPACE myts2 OWNER ics LOCATION '/home/space2';



-- ### 4. Создать БД
-- 1. Создать БД mydb
DROP DATABASE IF EXISTS mydb;
CREATE DATABASE mydb
    WITH
    OWNER = ics
    ENCODING = 'UTF8'
    TABLESPACE = myts1;

-- 2. Создать БД mytest
DROP DATABASE IF EXISTS mytest;
CREATE DATABASE mytest
    WITH
    OWNER = ics
    ENCODING = 'UTF8'
    TABLESPACE = myts2;

-- 3. Уничтожить БД mytest
DROP DATABASE IF EXISTS mytest;



-- ### 5. Создать схемы
-- - 1. myschem1 для mydb1
CREATE SCHEMA IF NOT EXISTS myschem1;

-- - 1. myschem2 для mydb1
CREATE SCHEMA IF NOT EXISTS myschem2;

-- - 3. Определитб текущую схему
select current_schema();
-- или
SHOW search_path;

-- - 4. Сделать myschem2 текущей
SET search_path TO myschem2;
SHOW search_path;



-- ### 6. Создать последовательность
-- переключаемся на схему 1
SET search_path TO myschem1;

-- 1. MySq c шагом 1
CREATE SEQUENCE IF NOT EXISTS mysq1
    INCREMENT 1
    START 1;


-- ### 7. Создать Тип
-- 1.Создать новый тип с описанием работника
CREATE TYPE fio AS
(
    name    character(40),
    soname  character(40),
    family  character(40),
    gender  character(1)
);



-- ### 8. Создать домен
/*
In PostgreSQL, a domain is a data type with optional constraints e.g.,
NOT NULL, CHECK etc. A domain has a unique name within the schema scope.
Domains are useful for centralizing management of fields with the common constraints.
*/
CREATE DOMAIN mydom AS
    INT CHECK (value < 7);



-- ### 9. Добавить расширение pg_freespacemap
CREATE EXTENSION pg_freespacemap;



-- ### 10.Таблицы
-- - a. Создать в схеме myschem1 таблицы согласно перечисленным ниже требованиям
-- - b. Создать связи между таблицами, как показано на логической схеме БД
-- - c. Посмотреть каталог где располагается таблица

/*
-----------------------------
| Заказ с аптечного склада  |  <- tit_in aka title_incoming
-----------------------------
|             | nom_in 1    |  <- полное описание единицы товара - в nomenclatura
|спецификация | nom_in 2    |
|             | nom_in 1    |
-----------------------------
*/
drop view if exists spec_report;
ALTER SEQUENCE mysq1 RESTART WITH 1;

DROP TABLE IF EXISTS head;
DROP TABLE IF EXISTS nom_in;
DROP TABLE IF EXISTS nom_out;
DROP TABLE IF EXISTS tit_out;
DROP TABLE IF EXISTS tit_in;
DROP TABLE IF EXISTS employers;
DROP TABLE IF EXISTS nomenclatura;
DROP TABLE IF EXISTS instruction;


/*
https://akela.mendelu.cz/~xcervena/IS/is/PowerDesigner-6/Definition%20files/sybase.def
http://infocenter-archive.sybase.com/help/index.jsp?topic=/com.sybase.stf.powerdesigner.docs_12.0.0/html/cdug/cdugp221.htm
I       = int
N%n     = numeric
SI      = smallint
N%s,%p  = float
A%n     = char(%n)
LVA     = varchar(255)
VA%n    = varchar(%n)
D       = datetime
TS      = timestamp
*/


-- - Таблица EMPLOYERS
--     - 1. Первичные ключи объявить на основе последовательности.
--     - 2. aSTelephone объявить как массив строк.
--     - 3. Добавить столбец FIO1 типа fio
CREATE TABLE employers
(
    ID_EMP              INT PRIMARY KEY DEFAULT nextval('mysq1'),
    SName               varchar(40),
    SFamily             varchar(40),
    FIO1                fio,
    SPosition           varchar(40),
    SSex                int,
    aSTelephone         text [],
    SOrganization       varchar(255),
    S_FIO_EMPL          varchar(50),
    SLogin              varchar(40),
    SPSW                char(10),
    IAccess             int,
    sTCHF               varchar(20) -- tech field
);

-- - Таблица NOMECLATURA
--     1. Первичные ключи объявить как автоинкрементальный (Serial)..
--     2. SCODE unique
--     3. dtINPOUT текущая дата  // такого вообще нет. поговорили с ним. решили убрать атрибут
CREATE TABLE nomenclatura
(
    ID_DRG               SERIAL PRIMARY KEY,
    REM_ID               numeric(6),
    sCODE                varchar(20) UNIQUE,
    DRUG_BARF            numeric(14),
    DRUG_NAME            varchar(255),
    PACK1_QTTY           numeric(4),
    NOM_QTTY             numeric(4),
    INVAL_DATE           date,
    siACTUAL             SMALLINT,
    sINSTRUCTION         varchar(255),
    CHECK_DATE           date,
    sNOTE                varchar(255),
    sTCHF                varchar(20)        -- tech field
);
comment on table nomenclatura is 'общий справочник товара';

-- - Таблица требований TIT_OUT
--     - 1. Первичные ключи объявить на основе автоинкрементальный (Serial).
--     - 2. iCODE уникальное поле
CREATE TABLE tit_out
(
    ID                   SERIAL PRIMARY KEY,
    iDIVISION            int,
    iCODE                numeric(6) UNIQUE,
    dDATE                date,      -- срок годности
    dLAST_CORR           date,
    iWORK_UP             int,
    dWORK_UP             timestamp WITHOUT TIME ZONE,
    sSUBSCR_APT          varchar(30),
    sSUBSCR_DIV          varchar(30),
    sNOTE                varchar(255),
    sTCHF                varchar(200)    -- tech field. в 20 символов не помещается ничего -> 200
);
comment on table tit_out is 'Расход';

-- - Таблица TIT_IN
--     1. Первичные ключи объявить автоинкрементальный
--     2. iCODE уникальный
CREATE TABLE tit_in
(
    ID_TIT               SERIAL PRIMARY KEY,
    EMP_ID               int,
    iCODE                varchar(16) UNIQUE,
    dDATE                date,      -- срок годности
    dLASTCORR            date,      -- последняя корректировка
    dWORKUP              timestamp WITHOUT TIME ZONE,
    nSUM                 numeric(11, 2),
    nCORR                numeric(11, 2),
    sSUBSCR              varchar(30),
    sNOTE                varchar(255),
    sTCHF                varchar(200),

    CONSTRAINT fk_emp
        FOREIGN KEY (EMP_ID)
        REFERENCES employers (ID_EMP)
);
comment on table tit_in is 'Приход';

-- - Таблицы NOM_IN
--     1. Первичные ключи определиться как автоинкрементальный (Serial).
--     2. Предусмотреть поля для внешних ключей ID_TIT и ID_NOM
--     3. Дата текущая
--     4. Связать таблицы как указано на схеме
CREATE TABLE nom_in
(
    ID_CL                SERIAL PRIMARY KEY,
    nQUANTITY            numeric(12, 3),
    nQUANT               numeric(12, 3),
    dLASTCORR            DATE NOT NULL DEFAULT CURRENT_DATE,
    dDATE                DATE NOT NULL DEFAULT CURRENT_DATE,    -- срок годности
    PACK1_QTTY           numeric(4),
    NOM_QTTY             numeric(4),
    siERROR              SMALLINT,
    sTCHF                varchar(200),
    ID_TIT               INT,   --nullable
    ID_NOM               INT,   --nullable

    CONSTRAINT fk_tit_in
        FOREIGN KEY (ID_TIT)
        REFERENCES tit_in (ID_TIT),
    CONSTRAINT fk_nomenclatura
        FOREIGN KEY (ID_NOM)
        REFERENCES nomenclatura (ID_DRG)
);
comment on table nom_in is 'накладная входящая. ака спецификация';

-- - Таблица NOM_OUT
--     1. Первичные ключи определиться как автоинкрементальный (Serial)
--     2. Дата по умолчанию текущая
--     3. Предусмотреть поля для внешних ключей ID_TIT и ID_NOM
--     4. Связать таблицы как указано на схеме
CREATE TABLE nom_out
(
    ID                   SERIAL PRIMARY KEY,
    nQUANT_C             numeric(12, 3),
    nQUANT_F             numeric(12, 3),
    dLAST_COOR           DATE NOT NULL DEFAULT CURRENT_DATE,
    siERROR              SMALLINT,
    dDATE                DATE NOT NULL DEFAULT CURRENT_DATE,
    sUNIT                char(1),
    nPRICE               numeric(12, 2),
    siPROPIS             SMALLINT,
    TIT_ID               numeric(10),
    sTCHF                varchar(200),
    ID_TIT               INT,    --nullable
    ID_NOM               INT,    --nullable

    CONSTRAINT fk_tit_out
        FOREIGN KEY (ID_TIT)
        REFERENCES tit_out (ID),
    CONSTRAINT fk_nomenclatura
        FOREIGN KEY (ID_NOM)
        REFERENCES nomenclatura (ID_DRG)
);
comment on table nom_out is 'накладная исходящая. aka специфкация';


-- - Таблица HEAD
--     1. Наследованием EMPLOYES
--     2. Дополнительное поле nadd - надбавка
CREATE TABLE head
(
	nadd                INT,
	PRIMARY KEY (ID_EMP)
)
INHERITS (employers);

-- - Таблица INSTRUCTION
--     1. поле для текстовой информации inst
--     2. поле ins поле типа tsvector
--
-- the tsvector data type, where ts stands for "text search");
-- to_tsquery for querying the vector for occurrences of certain words or phrases.
CREATE TABLE instruction
(
	inst                text,
	ins                 tsvector
);


-- - Временная таблица Temp EMPLOYERS
--     Структура таблицы такая же как и EMPLOYERS
--
-- Временная таблица существует только в рамках текущей сессии.
-- Если создать новое подключение и выплнить `SELECT * FROM temp_employers;`, то таблица не найдется
DROP TABLE IF EXISTS temp_employers;
CREATE TEMPORARY TABLE temp_employers
(
    ID_EMP              BIGINT PRIMARY KEY DEFAULT nextval('mysq1'),
    SName               varchar(40),
    SFamily             varchar(40),
    FIO1                fio,
    SPosition           varchar(40),
    SSex                int,
    aSTelephone         text [],
    SOrganization       varchar(255),
    S_FIO_EMPL          varchar(50),
    SLogin              varchar(40),
    SPSW                char(10),
    IAccess             int,
    sTCHF               varchar(200) -- tech field
);


-- Создать индекс по коду накладной, и дате поступления.
CREATE INDEX idx_nomanclatura_scode_inval_date
    ON nomenclatura(sCODE, INVAL_DATE);




-- ### 11. SQL-операторы
-- Insert

-- 1. Добавьте двух работников в EMPLOYERS; один c двумя телефонами,
-- один c тремя телефонами, заполнитть столбец FIO1
INSERT INTO employers(SName, SPosition, FIO1, aSTelephone) VALUES
('steve', 'worker', ('steve', 'rogers', 'marvel', 'm'), '{"8-800-555-3535", "8-915-144-4838"}'),
('peter', 'worker', ('peter', 'quill', 'guardians', 'm'), '{"8-123-358-4387"}'),
('tony', 'manager', ('tony', 'stark', 'marvel', 'm'), '{"+7-128-547-5491", "8-912-478-1928", "8-944-123-8549"}');

--2. Добавьте две номенклатуры в NOMECLATURA
INSERT INTO nomenclatura(sCODE, DRUG_NAME) VALUES
('12', 'арбидол'),
('23', 'смекта'),
('34', 'нурофен');


--3. Добавьте накладную без спецификации.
--  накладная без спецификации имеет только "заказ", но не имеет полей с единицами товара
-- delete from tit_in;
insert into tit_in(EMP_ID, sNOTE, dWORKUP) values
(1, 'Заказ для поликлиники #1488 г. Москва', CURRENT_TIMESTAMP),
(1, 'Заказ поликлиники МГТУ им бэтмана', CURRENT_TIMESTAMP),
(1, 'Пустой заказ без спецификации', CURRENT_TIMESTAMP);


--4. Добавьте номенклатуру в спецификацию.
-- добваляем единицы номенклатуры (как спецификацию) в накладную выше
-- delete from nom_in;
insert into nom_in(nQUANTITY, ID_TIT, ID_NOM, sTCHF, nQUANT) values
(20, 1, 1, '20 арбидолов для поликлиники 1488', 399),
(300, 1, 2, '300 смект по 89 руб, в пк 1488 все очень плохо', 89),
(500, 2, 3, 'дешевые нурофены для лучшего технического во вселенной', 39);


--5. Добавьте текст инструкции для обоих записей
-- инструкции есть только в номенклатуре
update nomenclatura
set sINSTRUCTION = '1 таблетку 2 раза в день'
where ID_DRG in (1, 3);
update nomenclatura
set sINSTRUCTION = 'размешать с водой, выпить'
where ID_DRG in (2);


--6. Создать функцию для дублирования строки с заданной номенклатурой.
CREATE OR REPLACE FUNCTION duplicate_nomenclatura(record_id int)
  RETURNS void AS
$func$
BEGIN
   EXECUTE format('INSERT INTO nomenclatura (
                     DRUG_NAME, sINSTRUCTION
                   )
                   SELECT DRUG_NAME, sINSTRUCTION FROM nomenclatura
                   WHERE ID_DRG = %s', record_id);
END
$func$ LANGUAGE plpgsql;
-- Select duplicate_nomenclatura(2);

--1. Сформировать список всех работников.
select SName, SPosition
from employers
where sPosition = 'worker';

--2. Выдать список всех руководителей
select SName, SPosition
from employers
where sPosition = 'manager';



-- #Отчеты по накладным
--1. Написать оператор Select который формирует список работников с данным телефоном первым в списке.
select SName, SPosition, aSTelephone
from employers
where 1 = array_position(aSTelephone, '8-800-555-3535');

--2. Написать оператор Select который формирует отчет обо всех накладных и спецификации,
-- относящиеся к ним. (Шапка должна включать номер товарной накладной, дату формирования
-- и название номенклатуры, количество, цену.) Создать Vew
create view spec_report as
select t.ID_TIT, t.sNOTE, t.dWORKUP, n.DRUG_NAME, ni.nQUANTITY, ni.nQUANT
from tit_in t
left join nom_in ni on t.ID_TIT = ni.ID_TIT
left join nomenclatura n on ni.ID_NOM = n.ID_DRG;

--3. Написать оператор Select который формирует накладные, имеющие спецификаци_и.
select t.ID_TIT, t.sNOTE
from tit_in t
where t.ID_TIT in (select ni.ID_TIT from nom_in ni);

--4. Написать оператор Select который формирует накладные, имеющих спецификаци_ю.
select t.ID_TIT, t.sNOTE
from tit_in t
where t.ID_TIT not in (select ni.ID_TIT from nom_in ni);

--5. Написать оператор Select который формирует спецификации.
select
--6. Написать оператор Select который формирует препараты, имеющие максимальную стоимость.
--7. Написать оператор Select который формирует накладными .
--отчет обо всех накладных, не отчет обо всех накладных, отчет с суммарной стоимостью
--отчет со списком накладных, отчет с товарными
--отчет о средней стоимости
--8. Написать оператор Select который формирует номенклатуры накладных на внутреннее перемещение.
--9. Написать оператор Select который формирует отчет о средней стоимости номенклатуры внутри каждой товарной накладной.
--10. Написать оператор Select который формирует отчет с перечнем накладных на имеющих минимальную и минимальную стоимость спецификации


### 12. Посмотреть план запросов
Для Select из п.11 первое задание


### 13.Создать функцию

- Функцию, которая случайным образом создает требования

что такое требования????
- Создать функцию, которая увеличивает стоимость номенклатуры на 30 процентов.
