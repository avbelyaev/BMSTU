
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
('34', 'нурофен'),
('56', 'сеалекс');


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
-- delete from nom_in;
insert into nom_in(nQUANTITY, ID_TIT, ID_NOM, sTCHF, nQUANT) values
(20, 1, 1, 'арбидолы для поликлиники', 399),
(300, 1, 2, 'смекты по 89 руб', 89),
(500, 2, 3, 'дешевые нурофены для лучшего технического', 39),
(20, 4, 4, 'для вмк', 299);

--5. Добавьте текст инструкции для обоих записей
-- инструкции есть только в номенклатуре
update nomenclatura
set sINSTRUCTION = '1 таблетку 2 раза в день'
where ID_DRG in (1, 3, 4);
update nomenclatura
set sINSTRUCTION = 'размешать с водой, выпить'
where ID_DRG in (2);

select * from nomenclatura;

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
Select duplicate_nomenclatura(2);

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
create view spec_overall_report as
select t.ID_TIT, t.sNOTE, t.dWORKUP, n.DRUG_NAME, ni.nQUANTITY, ni.nQUANT
from tit_in t
left join nom_in ni on t.ID_TIT = ni.ID_TIT
left join nomenclatura n on ni.ID_NOM = n.ID_DRG;

--3. Написать оператор Select который формирует отчет обо всех накладных, не имеющих спецификациии.
select t.ID_TIT, t.sNOTE
from tit_in t
where t.ID_TIT not in (select ni.ID_TIT from nom_in ni);

--4. Написать оператор Select который формирует отчет обо всех накладных, имеющих спецификацию.
select t.ID_TIT, t.sNOTE
from tit_in t
where t.ID_TIT in (select ni.ID_TIT from nom_in ni);

--5. Написать оператор Select который формирует отчет с суммарной стоимостью спецификации.
-- формируем спецификации для заказа 2
create view spec_prices as
select t.sNOTE, t.ID_TIT, sum(ni.nQUANTITY * ni.nQUANT) as price
from nom_in ni
join tit_in t on ni.ID_TIT = t.ID_TIT
group by t.sNOTE, t.ID_TIT;

--6. Написать оператор Select который формирует отчет со списком накладных, имеющих максимальную стоимость
select * from spec_prices
order by price DESC
limit 1;

--7. Написать оператор Select который формирует отчет с товарными накладными
select * from  nom_in;

--8. Написать оператор Select который формирует отчет о средней стоимости номенклатуры накладных на внутреннее перемещение.
select avg(price) from spec_prices;

--9. Написать оператор Select который формирует отчет о средней стоимости номенклатуры внутри каждой товарной накладной.
select t.sNOTE, avg(ni.nQUANTITY * ni.nQUANT) as price
from nom_in ni
join tit_in t on ni.ID_TIT = t.ID_TIT
group by t.sNOTE;

--10. Написать оператор Select который формирует отчет с перечнем накладных имеющих максимальную
select t.sNOTE, sum(ni.nQUANTITY * ni.nQUANT) as price
from nom_in ni
join tit_in t on ni.ID_TIT = t.ID_TIT
group by t.sNOTE
order by price DESC
limit 1;
-- и минимальную стоимость спецификации
select t.sNOTE, sum(ni.nQUANTITY * ni.nQUANT) as price
from nom_in ni
join tit_in t on ni.ID_TIT = t.ID_TIT
group by t.sNOTE
order by price ASC
limit 1;


-- ### Updata (Изменения значения строк)
-- 1. Изменить отчество в столбец FIO1 в первой записи
select * from employers;

update employers
set fio1.family = 'wozniak'
where ID_EMP = 1;

-- 2. Увеличить цену каждой номенклатуры на десять процентов.
select * from nom_in;

update nom_in
set nQUANT = nQUANT * 1.1;

-- 3. Переместить одну номенклатурную единицу из одной товарной накладной в другую.
-- Товарные накладные выбрать по своему усмотрению.
-- перемещение товара между накладными есть смена ID_TIT товара. в качестве товара берем 1ю единицу
update nom_in
set ID_TIT = 2
where ID_CL = 1;

-- 4. Поменять спецификации двух произвольных накладных на внутреннее перемещение.
select * from nom_in;

update nom_in
set ID_TIT = 2
where ID_TIT = 1;

update nom_in
set ID_TIT = 1
where ID_TIT = 2;


-- ### Удаление
-- 1. Удалить товарные накладные с минимальной и максимальной стоимостью спецификации.
-- === пересоздаем данные ===
select  * from spec_prices;

delete from tit_in
where ID_TIT in (
    select ID_TIT
    from spec_prices
    where price in (
        (select max(price) from spec_prices),
        (select min(price) from spec_prices)
    )
);


-- восстанавливаем данные (привязываем номенклатуу к другому заказу)
update nom_in
set ID_TIT = 1
WHERE ID_TIT is null;


-- 2. Удалить товарные накладные на сумму большую, чем средняя стоимость.
delete from tit_in
where ID_TIT in (
    select ID_TIT
    from spec_prices
    where price > (select avg(price) from spec_prices)
);

-- 3. Удалить товарные накладные без спецификации.
delete from tit_in
where ID_TIT not in (
    select ID_TIT
    from nom_in
);

-- 4. Удалить товарные накладные со спецификацией.
delete from tit_in
where ID_TIT in (
    select ID_TIT
    from nom_in
);



-- ### 12. Посмотреть план запросов
-- Для Select из п.11 первое задание
-- TODO explain



-- ### 13.Создать функцию
-- 1. Функцию, которая случайным образом создает требования (aka tit_in).
select * from tit_in;

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

-- select random_tit_in();


-- 2. Создать функцию, которая увеличивает стоимость номенклатуры на 30 процентов.
select * from nom_in;

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

select increase_nom_in_price_by_30_percent(4);


-- ### 14.Создать триггер
-- PostgreSQL only allows the execution of a user-defined function for the triggered action.
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
select * from tit_out;
select * from nom_out;

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

delete from tit_out
where id = 2;



-- ### 15.Создать Правило
-- Before the query is optimized, a rule can either replace the query with a different one or add additional queries.
-- These are then planned and executed instead of or together with the original query.

select * from nom_out;

-- Который реагирует на оператор insert с TIT_OUT и удаляющие все проасоциированые с этой строкой строки одной транзакцией.
-- TODO wtf is that
-- create rule rule_remove_associated
-- as on delete to tit_out
-- do instead blablabla



-- ### 16. FTS
-- Таблица INSTRUCTION
-- 1. C полем для текстовой информации inst и с поле ins поле типа tsvector
-- 2. Построить, поисковые вектора для всех строк таблицы Instruction
-- 3. Построить веса

-- the tsvector data type, where ts stands for "text search");
-- to_tsquery for querying the vector for occurrences of certain words or phrases.
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

select * from instruction;

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

-- запрос 'batman & joker' - видим как меняется вес в зависимости от положения бэтмена в строке
select inst, ins, ts_rank(ins, to_tsquery('joker & batman')) as rank
from instruction
order by rank desc;

-- запрос 'batman | joker' - наличие 2х слабых вариантов (D+C) перебило 1 сильный (B)
select inst, ins, ts_rank(ins, to_tsquery('joker | batman')) as rank
from instruction
order by rank desc;
