
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
SET search_path TO myschem1;

drop view if exists spec_overall_report;
drop view if exists spec_prices;
ALTER SEQUENCE mysq1 RESTART WITH 1;

DROP TABLE IF EXISTS head;
DROP TABLE IF EXISTS nom_in;
DROP TABLE IF EXISTS nom_out;
DROP TABLE IF EXISTS tit_out;
DROP TABLE IF EXISTS tit_in;
DROP TABLE IF EXISTS employers;
DROP TABLE IF EXISTS nomenclatura;


/*
Конвертация Sybase'ных типов в нормальные типы данных

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
--     3. dtINPOUT текущая дата  // такого вообще нет. решили убрать атрибут
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
comment on table tit_out is 'Расход. ака накладная';

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
comment on table tit_in is 'Приход. ака накладная';

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
        REFERENCES tit_in (ID_TIT)
        ON DELETE set null,
    CONSTRAINT fk_nomenclatura
        FOREIGN KEY (ID_NOM)
        REFERENCES nomenclatura (ID_DRG)
);
comment on table nom_in is 'спецификация входящей накладной';

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
comment on table nom_out is 'специфкация исходящей накладной';


-- - Таблица HEAD
--     1. Наследованием EMPLOYES
--     2. Дополнительное поле nadd - надбавка
CREATE TABLE head
(
	nadd                INT,
	PRIMARY KEY (ID_EMP)
)
INHERITS (employers);


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

-- проверям, что все таблицы созданы
select table_catalog, table_schema, table_name
from information_schema.tables
where table_schema = 'myschem1';
--  table_catalog | table_schema |  table_name
-- ---------------+--------------+--------------
--  mydb          | myschem1     | nomenclatura
--  mydb          | myschem1     | nom_in
--  mydb          | myschem1     | tit_in
--  mydb          | myschem1     | tit_out
--  mydb          | myschem1     | employers
--  mydb          | myschem1     | nom_out
--  mydb          | myschem1     | head
-- (7 rows)

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
('34', 'алко-зельтцер'),
('56', 'андрогель');


--3. Добавьте накладную без спецификации.
--  накладная без спецификации имеет только "заказ", но не имеет полей с единицами товара
-- delete from tit_in;
insert into tit_in(EMP_ID, sNOTE, dWORKUP) values
(1, 'Заказ для поликлиники #1488 г. Москва', CURRENT_TIMESTAMP),
(1, 'Заказ поликлиники МГТУ им бэтмана', CURRENT_TIMESTAMP),
(1, 'Пустой заказ без спецификации', CURRENT_TIMESTAMP),
(1, 'Заказ для МГУ', CURRENT_TIMESTAMP);

--4. Добавьте номенклатуру в спецификацию.
-- добваляем единицы номенклатуры (как спецификацию) в накладную выше

/*
 ввиду сложности и непривычности модели данных, далее сделано следующее упрощение:
 - под nQUANTITY будет подразумевать КОЛИЧЕСВТО ТОВАРА
 - под nQUANT - ЦЕНУ ЗА ЕДИНИЦУ ОТВАРА
 т.к. необходимо показать именно работу с базой данных, а не знания в области логистики
 */
insert into nom_in(nQUANTITY, ID_TIT, ID_NOM, sTCHF, nQUANT) values
(20, 1, 1, 'арбидолы для поликлиники', 399),
(300, 1, 2, 'смекты по 89 руб', 89),
(500, 2, 3, 'препараты для лучшего технического', 39),
(20, 4, 4, 'все лучшее для вмк', 299);

--5. Добавьте текст инструкции для обоих записей
-- инструкции есть только в номенклатуре
update nomenclatura
set sINSTRUCTION = '1 таблетку 2 раза в день'
where ID_DRG in (1, 3, 4);
update nomenclatura
set sINSTRUCTION = 'размешать с водой, выпить'
where ID_DRG in (2);

select id_drg, scode, drug_barf, drug_name, sinstruction, stchf from nomenclatura;
--  id_drg | scode | drug_barf |   drug_name   |       sinstruction        | stchf
-- --------+-------+-----------+---------------+---------------------------+-------
--       1 | 12    |           | арбидол       | 1 таблетку 2 раза в день  |
--       3 | 34    |           | алко-зельтцер | 1 таблетку 2 раза в день  |
--       4 | 56    |           | андрогель     | 1 таблетку 2 раза в день  |
--       2 | 23    |           | смекта        | размешать с водой, выпить |
-- (4 rows)

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

-- call
Select duplicate_nomenclatura(2);

-- проверка
select id_drg, scode, drug_barf, drug_name, sinstruction, stchf
from nomenclatura;
--  id_drg | scode | drug_barf |   drug_name   |       sinstruction        | stchf
-- --------+-------+-----------+---------------+---------------------------+-------
--       1 | 12    |           | арбидол       | 1 таблетку 2 раза в день  |
--       3 | 34    |           | алко-зельтцер | 1 таблетку 2 раза в день  |
--       4 | 56    |           | андрогель     | 1 таблетку 2 раза в день  |
--       2 | 23    |           | смекта        | размешать с водой, выпить |
--       5 |       |           | смекта        | размешать с водой, выпить |
-- (5 rows)

--1. Сформировать список всех работников.
select SName, SPosition
from employers
where sPosition = 'worker';
--  sname | sposition
-- -------+-----------
--  steve | worker
--  peter | worker
-- (2 rows)

--2. Выдать список всех руководителей
select SName, SPosition
from employers
where sPosition = 'manager';
--  sname | sposition
-- -------+-----------
--  tony  | manager
-- (1 row)



-- #Отчеты по накладным
--1. Написать оператор Select который формирует список работников с данным телефоном первым в списке.
select SName, SPosition, aSTelephone
from employers
where 1 = array_position(aSTelephone, '8-800-555-3535');
--  sname | sposition |           astelephone
-- -------+-----------+---------------------------------
--  steve | worker    | {8-800-555-3535,8-915-144-4838}
-- (1 row)

--2. Написать оператор Select который формирует отчет обо всех накладных и спецификации,
-- относящиеся к ним. (Шапка должна включать
-- - номер товарной накладной
-- - (добавляем так же sNOte для наглядности данных)
-- - дату формирования
-- - название номенклатуры
-- - количество
-- - цену.
-- Создать Vew
create view spec_overall_report as
select t.ID_TIT, t.sNOTE, t.dWORKUP, n.DRUG_NAME, ni.nQUANTITY, ni.nQUANT
from tit_in t
left join nom_in ni on t.ID_TIT = ni.ID_TIT
left join nomenclatura n on ni.ID_NOM = n.ID_DRG;

select * from spec_overall_report;
--  id_tit |                 snote                 |          dworkup          |   drug_name   | nquantity | nquant
-- --------+---------------------------------------+---------------------------+---------------+-----------+---------
--       1 | Заказ для поликлиники #1488 г. Москва | 2019-10-28 09:50:29.82451 | арбидол       |    20.000 | 399.000
--       1 | Заказ для поликлиники #1488 г. Москва | 2019-10-28 09:50:29.82451 | смекта        |   300.000 |  89.000
--       2 | Заказ поликлиники МГТУ им бэтмана     | 2019-10-28 09:50:29.82451 | алко-зельтцер |   500.000 |  39.000
--       4 | Заказ для МГУ                         | 2019-10-28 09:50:29.82451 | андрогель     |    20.000 | 299.000
--       3 | Пустой заказ без спецификации         | 2019-10-28 09:50:29.82451 |               |           |
-- (5 rows)

--3. Написать оператор Select который формирует отчет обо всех накладных, не имеющих спецификациии.
select t.ID_TIT, t.sNOTE
from tit_in t
where t.ID_TIT not in (select ni.ID_TIT from nom_in ni);
--  id_tit |             snote
-- --------+-------------------------------
--       3 | Пустой заказ без спецификации
-- (1 row)

--4. Написать оператор Select который формирует отчет обо всех накладных, имеющих спецификацию.
select t.ID_TIT, t.sNOTE
from tit_in t
where t.ID_TIT in (select ni.ID_TIT from nom_in ni);
--  id_tit |                 snote
-- --------+---------------------------------------
--       1 | Заказ для поликлиники #1488 г. Москва
--       2 | Заказ поликлиники МГТУ им бэтмана
--       4 | Заказ для МГУ
-- (3 rows)

--5. Написать оператор Select который формирует отчет с суммарной стоимостью спецификации.
-- формируем спецификации для заказа 2
create view spec_prices as select
       t.sNOTE,
       t.ID_TIT,
       sum(coalesce(ni.nQUANTITY * ni.nQUANT, 0)) as overall_price,
       avg(coalesce(ni.nQUANTITY * ni.nQUANT, 0)) as avg_price
from nom_in ni
right join tit_in t on ni.ID_TIT = t.ID_TIT
group by t.sNOTE, t.ID_TIT;

select * from spec_prices;
--                  snote                 | id_tit | overall_price |       avg_price
-- ---------------------------------------+--------+---------------+------------------------
--  Заказ для МГУ                         |      4 |   5980.000000 |  5980.0000000000000000
--  Заказ поликлиники МГТУ им бэтмана     |      2 |  19500.000000 | 19500.0000000000000000
--  Пустой заказ без спецификации         |      3 |             0 | 0.00000000000000000000
--  Заказ для поликлиники #1488 г. Москва |      1 |  34680.000000 |     17340.000000000000
-- (4 rows)

--6. Написать оператор Select который формирует отчет со списком накладных, имеющих максимальную стоимость
select * from spec_prices
order by overall_price DESC
limit 1;
--                  snote                 | id_tit | overall_price |     avg_price
-- ---------------------------------------+--------+---------------+--------------------
--  Заказ для поликлиники #1488 г. Москва |      1 |  34680.000000 | 17340.000000000000
-- (1 row)

--7. Написать оператор Select который формирует отчет с товарными накладными
select id_cl, nquantity, nquant, ddate, stchf, id_tit, id_nom
from  nom_in;
--  id_cl | nquantity | nquant  |   ddate    |               stchf                | id_tit | id_nom
-- -------+-----------+---------+------------+------------------------------------+--------+--------
--      1 |    20.000 | 399.000 | 2019-10-28 | арбидолы для поликлиники           |      1 |      1
--      2 |   300.000 |  89.000 | 2019-10-28 | смекты по 89 руб                   |      1 |      2
--      3 |   500.000 |  39.000 | 2019-10-28 | препараты для лучшего технического |      2 |      3
--      4 |    20.000 | 299.000 | 2019-10-28 | все лучшее для вмк                 |      4 |      4
-- (4 rows)

--8. Написать оператор Select который формирует отчет о средней стоимости номенклатуры накладных на внутреннее перемещение.
select avg_price from spec_prices;
--        avg_price
-- ------------------------
--   5980.0000000000000000
--  19500.0000000000000000
--  0.00000000000000000000
--      17340.000000000000
-- (4 rows)

--9. Написать оператор Select который формирует отчет о средней стоимости номенклатуры внутри каждой товарной накладной.
select t.sNOTE, avg(coalesce(ni.nQUANTITY * ni.nQUANT, 0)) as price
from nom_in ni
right join tit_in t on ni.ID_TIT = t.ID_TIT
group by t.sNOTE;
--                  snote                 |         price
-- ---------------------------------------+------------------------
--  Пустой заказ без спецификации         | 0.00000000000000000000
--  Заказ поликлиники МГТУ им бэтмана     | 19500.0000000000000000
--  Заказ для поликлиники #1488 г. Москва |     17340.000000000000
--  Заказ для МГУ                         |  5980.0000000000000000
-- (4 rows)

--10. Написать оператор Select который формирует отчет с перечнем накладных имеющих максимальную стоимость спецификации
select sNOTE, overall_price
from spec_prices
order by overall_price DESC
limit 1;
--                  snote                 |    price
-- ---------------------------------------+--------------
--  Заказ для поликлиники #1488 г. Москва | 34680.000000
-- (1 row)


-- и минимальную стоимость спецификации
select sNOTE, overall_price
from spec_prices
order by overall_price ASC
limit 1;
--              snote             | overall_price
-- -------------------------------+---------------
--  Пустой заказ без спецификации |             0
-- (1 row)


-- ### Update (Изменения значения строк)
-- 1. Изменить отчество в столбец FIO1 в первой записи
select ID_EMP, SName, FIO1, aSTelephone
from employers;
-- 1 | steve | ("steve ","rogers ","marvel     ",m) | {8-800-555-3535,8-915-144-4838}
-- 2 | peter | ("peter ","quill  ","guardians  ",m) | {8-123-358-4387}
-- 3 | tony  | ("tony  ","stark  ","marvel     ",m) | {+7-128-547-5491,8-912-478-1928,8-944-123-8549}
-- (3 rows)

update employers
set fio1.family = 'wozniak'
where ID_EMP = 1;
-- 2 | peter | ("peter  ","quill   ","guardians ",m)    | {8-123-358-4387}
-- 3 | tony  | ("tony  ","stark  ","marvel  ",m)        | {+7-128-547-5491,8-912-478-1928,8-944-123-8549}
-- 1 | steve | ("steve  ","rogers   ","wozniak ",m)     | {8-800-555-3535,8-915-144-4838}
-- (3 rows)

-- 2. Увеличить цену каждой номенклатуры на десять процентов.
-- смотрим изначальную цену
select id_cl, nquantity, nquant, stchf, id_tit, id_nom
from nom_in;
--  id_cl | nquantity | nquant  |               stchf                | id_tit | id_nom
-- -------+-----------+---------+------------------------------------+--------+--------
--      1 |    20.000 | 399.000 | арбидолы для поликлиники           |      1 |      1
--      2 |   300.000 |  89.000 | смекты по 89 руб                   |      1 |      2
--      3 |   500.000 |  39.000 | препараты для лучшего технического |      2 |      3
--      4 |    20.000 | 299.000 | все лучшее для вмк                 |      4 |      4
-- (4 rows)

-- увеличиваем
update nom_in
set nQUANT = nQUANT * 1.1;

-- смотрим новую цену
--  id_cl | nquantity | nquant  |               stchf                | id_tit | id_nom
-- -------+-----------+---------+------------------------------------+--------+--------
--      1 |    20.000 | 438.900 | арбидолы для поликлиники           |      1 |      1
--      2 |   300.000 |  97.900 | смекты по 89 руб                   |      1 |      2
--      3 |   500.000 |  42.900 | препараты для лучшего технического |      2 |      3
--      4 |    20.000 | 328.900 | все лучшее для вмк                 |      4 |      4
-- (4 rows)

-- 3. Переместить одну номенклатурную единицу из одной товарной накладной в другую.
-- Товарные накладные выбрать по своему усмотрению.
-- перемещение товара между накладными есть смена ID_TIT товара. в качестве товара берем 1ю единицу
update nom_in
set ID_TIT = 2
where ID_CL = 1;
--  id_cl | nquantity | nquant  |               stchf                | id_tit | id_nom
-- -------+-----------+---------+------------------------------------+--------+--------
--      2 |   300.000 |  97.900 | смекты по 89 руб                   |      1 |      2
--      3 |   500.000 |  42.900 | препараты для лучшего технического |      2 |      3
--      4 |    20.000 | 328.900 | все лучшее для вмк                 |      4 |      4
--      1 |    20.000 | 438.900 | арбидолы для поликлиники           |      2 |      1
-- (4 rows)

-- 4. Поменять спецификации двух произвольных накладных на внутреннее перемещение.
update nom_in
set ID_TIT = 2
where ID_TIT = 1;

update nom_in
set ID_TIT = 1
where ID_TIT = 2;
--  id_cl | nquantity | nquant  |               stchf                | id_tit | id_nom
-- -------+-----------+---------+------------------------------------+--------+--------
--      4 |    20.000 | 328.900 | все лучшее для вмк                 |      4 |      4
--      3 |   500.000 |  42.900 | препараты для лучшего технического |      1 |      3
--      1 |    20.000 | 438.900 | арбидолы для поликлиники           |      1 |      1
--      2 |   300.000 |  97.900 | смекты по 89 руб                   |      1 |      2
-- (4 rows)


-- ### Удаление
-- === ПЕРЕСОЗДАЕМ ДАННЫЕ для наглядности работы ===
select  * from spec_prices;
--                  snote                 | id_tit | overall_price |       avg_price
-- ---------------------------------------+--------+---------------+------------------------
--  Заказ для МГУ                         |      4 |   6578.000000 |  6578.0000000000000000
--  Заказ поликлиники МГТУ им бэтмана     |      2 |  21450.000000 |     21450.000000000000
--  Пустой заказ без спецификации         |      3 |             0 | 0.00000000000000000000
--  Заказ для поликлиники #1488 г. Москва |      1 |  38148.000000 |     19074.000000000000
-- (4 rows)

-- 1. Удалить товарные накладные с минимальной и максимальной стоимостью спецификации.
delete from tit_in
where ID_TIT in (
    select ID_TIT
    from spec_prices
    where overall_price in (
        (select max(overall_price) from spec_prices),
        (select min(overall_price) from spec_prices)
    )
);
--                snote               | id_tit | overall_price |       avg_price
-- -----------------------------------+--------+---------------+-----------------------
--  Заказ для МГУ                     |      4 |   6578.000000 | 6578.0000000000000000
--  Заказ поликлиники МГТУ им бэтмана |      2 |  21450.000000 |    21450.000000000000
-- (2 rows)

-- восстанавливаем данные (привязываем номенклатуу к другому заказу)
update nom_in
set ID_TIT = 4
WHERE ID_TIT is null;

select  * from spec_prices;
--                snote               | id_tit | overall_price |     avg_price
-- -----------------------------------+--------+---------------+--------------------
--  Заказ для МГУ                     |      4 |  44726.000000 | 14908.666666666667
--  Заказ поликлиники МГТУ им бэтмана |      2 |  21450.000000 | 21450.000000000000
-- (2 rows)

-- 2. Удалить товарные накладные на сумму большую, чем средняя стоимость.
-- здесь нужно именно среднее значение всех общих сумм
delete from tit_in
where ID_TIT in (
    select ID_TIT
    from spec_prices
    where overall_price > (select avg(overall_price) from spec_prices)
);
--                snote               | id_tit | overall_price |     avg_price
-- -----------------------------------+--------+---------------+--------------------
--  Заказ поликлиники МГТУ им бэтмана |      2 |  21450.000000 | 21450.000000000000
-- (1 row)

-- 3. Удалить товарные накладные без спецификации.
-- т.е. тайтлы (tit_in), ID которых нет среди входящей номенклатуры.
-- это, например, пустой заказ
delete from tit_in
where ID_TIT not in (
    select ID_TIT
    from nom_in
);

-- проверяем
select ID_TIT, sNOTE from tit_in;
--  id_tit |                 snote
-- --------+---------------------------------------
--       1 | Заказ для поликлиники #1488 г. Москва
--       2 | Заказ поликлиники МГТУ им бэтмана
--       4 | Заказ для МГУ
-- (3 rows)

-- 4. Удалить товарные накладные со спецификацией.
-- ожидаем, что удалятся все накладные, т.к. у всех наклданых есть спецификации
delete from tit_in
where ID_TIT in (
    select ID_TIT
    from nom_in
);

-- проверяем
select ID_TIT, sNOTE from tit_in;
--  id_tit | snote
-- --------+-------
-- (0 rows)



-- ### 12. Посмотреть план запросов
-- Для Select из п.11 первое задание

-- ПЕРЕСОЗДАЕМ ДАННЫЕ
-- берем следующий запрос в качесте примера:
select
       t.sNOTE,
       t.ID_TIT,
       sum(coalesce(ni.nQUANTITY * ni.nQUANT, 0)) as overall_price,
       avg(coalesce(ni.nQUANTITY * ni.nQUANT, 0)) as avg_price
from nom_in ni
right join tit_in t on ni.ID_TIT = t.ID_TIT
group by t.sNOTE, t.ID_TIT;
--  id_tit |                 snote                 |          dworkup           | nquantity | nquant
-- --------+---------------------------------------+----------------------------+-----------+---------
--       1 | Заказ для поликлиники #1488 г. Москва | 2019-10-28 11:31:33.733899 |    20.000 | 399.000
--       1 | Заказ для поликлиники #1488 г. Москва | 2019-10-28 11:31:33.733899 |   300.000 |  89.000
--       4 | Заказ для МГУ                         | 2019-10-28 11:31:33.733899 |    20.000 | 299.000
-- (3 rows)


explain (analyze) select t.ID_TIT, t.sNOTE, t.dWORKUP, ni.nQUANTITY, ni.nQUANT
from tit_in t
left join nom_in ni on t.ID_TIT = ni.ID_TIT
where ni.nQUANT > 50;

-- При выполнении запроса последовательно считывается каждая запись таблицы (Seq Scan)
/*
                           QUERY PLAN
----------------------------------------------------------------
 Nested Loop  (cost=0.00..2.14 rows=1 width=560)
   Join Filter: (t.id_tit = ni.id_tit)
   ->  Seq Scan on nom_in ni  (cost=0.00..1.05 rows=1 width=36)
         Filter: (nquant > '50'::numeric)
   ->  Seq Scan on tit_in t  (cost=0.00..1.04 rows=4 width=528)
(5 rows)
 */

-- принудительно включаем использование индекса, запретив Seq Scan
SET enable_seqscan TO off;
/*
                                         QUERY PLAN
-------------------------------------------------------------------------------------------
 Nested Loop  (cost=10000000000.13..10000000013.30 rows=1 width=560)
   Join Filter: (t.id_tit = ni.id_tit)
   ->  Index Scan using idx_tit_in_id_tit on tit_in t  (cost=0.13..12.19 rows=4 width=528)
   ->  Materialize  (cost=10000000000.00..10000000001.05 rows=1 width=36)
         ->  Seq Scan on nom_in ni  (cost=10000000000.00..10000000001.05 rows=1 width=36)
               Filter: (nquant > '50'::numeric)
(6 rows)
 */

-- как видим, стоимость запроса увеличилась
-- При выборке практически всей таблицы использование индекса
-- только увеличивает cost и время выполнения запроса. Планировщик не глуп.
-- возращаем seq scan на место
SET enable_seqscan TO off;


-- ### 13.Создать функцию

-- смотрим текщие наклданые
select t.ID_TIT, t.sNOTE, t.dWORKUP
from tit_in t;
--  id_tit |                 snote                 |          dworkup
-- --------+---------------------------------------+----------------------------
--       1 | Заказ для поликлиники #1488 г. Москва | 2019-10-28 11:31:33.733899
--       2 | Заказ поликлиники МГТУ им бэтмана     | 2019-10-28 11:31:33.733899
--       3 | Пустой заказ без спецификации         | 2019-10-28 11:31:33.733899
--       4 | Заказ для МГУ                         | 2019-10-28 11:31:33.733899
-- (4 rows)

-- 1. Функцию, которая случайным образом создает требования (aka tit_in).
CREATE OR REPLACE FUNCTION random_tit_in()
  RETURNS void AS
$func$
BEGIN
   EXECUTE format(
       'INSERT INTO tit_in (EMP_ID, sNOTE, dWORKUP) VALUES
        (%s, ''%s'', CURRENT_TIMESTAMP)',
       floor(random() * 3 + 1)::int,
       md5(random()::text)
       );
END
$func$ LANGUAGE plpgsql;

-- вызов
select random_tit_in();

-- смтрим что поменялось
select t.ID_TIT, t.sNOTE, t.dWORKUP
from tit_in t;
--  id_tit |                 snote                 |          dworkup
-- --------+---------------------------------------+----------------------------
--       1 | Заказ для поликлиники #1488 г. Москва | 2019-10-28 11:31:33.733899
--       2 | Заказ поликлиники МГТУ им бэтмана     | 2019-10-28 11:31:33.733899
--       3 | Пустой заказ без спецификации         | 2019-10-28 11:31:33.733899
--       4 | Заказ для МГУ                         | 2019-10-28 11:31:33.733899
--       5 | 86103e02284bb5c5f4563041ed6d8e84      | 2019-10-28 12:56:29.597853
-- (5 rows)


-- 2. Создать функцию, которая увеличивает стоимость номенклатуры на 30 процентов.
select id_cl, nquantity, nquant, stchf, id_tit, id_nom
from nom_in;
--  id_cl | nquantity | nquant  |               stchf                | id_tit | id_nom
-- -------+-----------+---------+------------------------------------+--------+--------
--      1 |    20.000 | 399.000 | арбидолы для поликлиники           |      1 |      1
--      2 |   300.000 |  89.000 | смекты по 89 руб                   |      1 |      2
--      3 |   500.000 |  39.000 | препараты для лучшего технического |      2 |      3
--      4 |    20.000 | 299.000 | все лучшее для вмк                 |      4 |      4
-- (4 rows)

CREATE OR REPLACE FUNCTION increase_nom_in_price_by_30_percent(id int)
  RETURNS void AS
$func$
BEGIN
   EXECUTE format(
       'UPDATE nom_in
        SET nQUANT = 1.3 * nQUANT
        WHERE id_cl = %s', id);
END
$func$ LANGUAGE plpgsql;

-- увеличиваем цену (nQuant) входной номенклатуры с ID = 4
select increase_nom_in_price_by_30_percent(4);
--  id_cl | nquantity | nquant  |               stchf                | id_tit | id_nom
-- -------+-----------+---------+------------------------------------+--------+--------
--      1 |    20.000 | 399.000 | арбидолы для поликлиники           |      1 |      1
--      2 |   300.000 |  89.000 | смекты по 89 руб                   |      1 |      2
--      3 |   500.000 |  39.000 | препараты для лучшего технического |      2 |      3
--      4 |    20.000 | 388.700 | все лучшее для вмк                 |      4 |      4
-- (4 rows)

-- ### 14.Создать триггер
-- PG разрешает лишь выполнение польовательских функций в качестве триггера действия
-- добавляем в tit_out/nom_out данные
insert into tit_out(dDATE, sSUBSCR_APT, dWORK_UP) values
(CURRENT_DATE, 'Выгрузка 228', CURRENT_TIMESTAMP),
(CURRENT_DATE, 'Выгрузка 322 пункта 1337', CURRENT_TIMESTAMP);

insert into nom_out(ID_TIT, ID_NOM, nPRICE, sTCHF) values
(1, 1, 350, 'арбидолы'),
(1, 2, 59, 'смекты'),
(2, 3, 240, 'нурофены');

-- Триггер реагирует на удаление строки из TIT_OUT.
-- Он должен удалять все проассоциированные с этой строкой строки одной транзакцией.
select id, ddate, dwork_up, ssubscr_apt
from tit_out;
--  id |   ddate    |          dwork_up          |       ssubscr_apt
-- ----+------------+----------------------------+--------------------------
--   1 | 2019-10-28 | 2019-10-28 12:00:11.540009 | Выгрузка 228
--   2 | 2019-10-28 | 2019-10-28 12:00:11.540009 | Выгрузка 322 пункта 1337
-- (2 rows)

select id, ddate, nprice, stchf, id_tit, id_nom
from nom_out;
--  id |   ddate    | nprice |  stchf   | id_tit | id_nom
-- ----+------------+--------+----------+--------+--------
--   1 | 2019-10-28 | 350.00 | арбидолы |      1 |      1
--   2 | 2019-10-28 |  59.00 | смекты   |      1 |      2
--   3 | 2019-10-28 | 240.00 | нурофены |      2 |      3
-- (3 rows)

-- в PG триггеимая функция по умолчанию вызвается в транзакции
CREATE OR REPLACE FUNCTION rm_nom_out() RETURNS trigger AS
$$BEGIN
    delete from nom_out
    where ID_TIT = OLD.id;
    return OLD;
END;$$ LANGUAGE plpgsql;

drop trigger if exists trigger_remove_associated on tit_out;

create trigger trigger_remove_associated
before delete on tit_out
for each row execute procedure rm_nom_out();

-- удаляем, вызывая триггер
delete from tit_out
where id = 2;

-- проверяем, что произошло удаление
select id, ddate, dwork_up, ssubscr_apt
from tit_out;
--  id |   ddate    |          dwork_up          | ssubscr_apt
-- ----+------------+----------------------------+--------------
--   1 | 2019-10-28 | 2019-10-28 12:00:11.540009 | Выгрузка 228
-- (1 row)

select id, ddate, nprice, stchf, id_tit, id_nom
from nom_out;
--  id |   ddate    | nprice |  stchf   | id_tit | id_nom
-- ----+------------+--------+----------+--------+--------
--   1 | 2019-10-28 | 350.00 | арбидолы |      1 |      1
--   2 | 2019-10-28 |  59.00 | смекты   |      1 |      2
-- (2 rows)



-- ### 15.Создать Правило
-- прежде, чем запрос будет оптимизирован, правило может изменить запрос на один или несколько других запросов
-- в итоге буду исполнены они, а не изнечальный запрос
-- можно сказать, что это макрос для запросов

select * from nom_out;

-- Который реагирует на оператор insert с TIT_OUT и удаляющие все проасоциированые с этой строкой строки одной транзакцией.
-- как это вообще понимать? условие некорректно

-- корректируем условия для наглядности:
drop rule if exists rule_select_instead_of_insert on tit_out;

create rule rule_select_instead_of_insert
as on delete to tit_out
do instead select 'fool ya! you cannot remove' from tit_out;

-- теперь при попытке удлить что-либо, правило макрос переписывать запрос на другой
delete from tit_out;
--           ?column?
-- ----------------------------
--  fool ya! you cannot remove
-- (1 row)



-- ### 16. FTS
-- Таблица INSTRUCTION
-- 1. C полем для текстовой информации inst и с поле ins поле типа tsvector

-- tsvector == "text search vector"
DROP TABLE IF EXISTS instruction;

CREATE TABLE instruction
(
	inst                text,
	ins                 tsvector
);

insert into instruction(inst, ins) values
('i1', to_tsvector(null)),
('i2', to_tsvector(null)),
('i3', to_tsvector(null)),
('i4', to_tsvector('Hello Darkness my old friend')),
('i5', to_tsvector('I''ve come to talk with you again'));

-- 2. Построить, поисковые вектора для всех строк таблицы Instruction
-- 3. Построить веса

-- в PG можно самом задать весь элемента вектора от A до D.
-- как правило, для названия статьи выбирают A. для ключевых слов - B. сама статья - C/D
-- максимальный вес джокера
UPDATE instruction
SET ins = setweight(to_tsvector('JOKER'), 'A')
where inst = 'i1';

-- большой вес джокера
UPDATE instruction
SET ins = setweight(to_tsvector('JoKeR'), 'B')
where inst = 'i2';

-- минимальный вес джокера, но есть бэтмен
UPDATE instruction SET ins =
    setweight(to_tsvector('jOKEr'), 'D')    ||
    setweight(to_tsvector('joker dark batman knight NOLAN bale'), 'C')
where inst = 'i3';

select * from instruction;
--  inst |                                 ins
-- ------+---------------------------------------------------------------------
--  i4   | 'dark':2 'friend':5 'hello':1 'old':4
--  i5   | 'come':3 'talk':5 've':2
--  i1   | 'joker':1A
--  i2   | 'joker':1B
--  i3   | 'bale':7C 'batman':4C 'dark':3C 'joker':1,2C 'knight':5C 'nolan':6C
-- (5 rows)

-- запрос 'batman & joker' - видим как меняется вес в зависимости от положения бэтмена в строке
select inst, ins, ts_rank(ins, to_tsquery('joker & batman')) as rank
from instruction
order by rank desc;
--  inst |                                 ins                                 |   rank
-- ------+---------------------------------------------------------------------+----------
--  i3   | 'bale':7C 'batman':4C 'dark':3C 'joker':1,2C 'knight':5C 'nolan':6C | 0.307563
--  i4   | 'dark':2 'friend':5 'hello':1 'old':4                               |    1e-20
--  i5   | 'come':3 'talk':5 've':2                                            |    1e-20
--  i1   | 'joker':1A                                                          |    1e-20
--  i2   | 'joker':1B                                                          |    1e-20
-- (5 rows)

-- запрос 'batman | joker' - наличие 2х слабых вариантов (D+C) перебило 1 более сильный (B)
select inst, ins, ts_rank(ins, to_tsquery('joker | batman')) as rank
from instruction
order by rank desc;
--  inst |                                 ins                                 |   rank
-- ------+---------------------------------------------------------------------+----------
--  i1   | 'joker':1A                                                          | 0.303964
--  i3   | 'bale':7C 'batman':4C 'dark':3C 'joker':1,2C 'knight':5C 'nolan':6C | 0.151982
--  i2   | 'joker':1B                                                          | 0.121585
--  i4   | 'dark':2 'friend':5 'hello':1 'old':4                               |        0
--  i5   | 'come':3 'talk':5 've':2                                            |        0
-- (5 rows)
