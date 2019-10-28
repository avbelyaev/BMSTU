-- noinspection SqlNoDataSourceInspectionForFile

-- ### 3. Создать tablespace
-- 1. Создать tablespace “myts1” (myts2) в директории MyDB1 (MyDB2) с пользователем postgres
/*
- в контейнере выпоняем:
mkdir -p /home/space1 && chown postgres /home/space1
mkdir -p /home/space2 && chown postgres /home/space2

- подчключаемся datagrip'ом, pgadmin'ом, либо PSQL'ем.
в случае psql:
psql --username ics --password ics --dbname ics
*/
-- - по отдельности запускаем:
create tablespace myts1 owner ics location '/home/space1';
create tablespace myts2 owner ics location '/home/space2';



-- ### 4. Создать БД
-- 1. Создать БД mydb
-- DROP DATABASE IF EXISTS mydb;
CREATE DATABASE mydb
    WITH
    OWNER = ics
    ENCODING = 'UTF8'
    TABLESPACE = myts1;

-- 2. Создать БД mytest
CREATE DATABASE mytest
    WITH
    OWNER = ics
    ENCODING = 'UTF8'
    TABLESPACE = myts2;

-- 3. Уничтожить БД mytest
DROP DATABASE IF EXISTS mytest;



-- ### 5. Создать схемы
/*
выбрать в pgadmin или datagrip'е конкретную БД с которой буем работать.
в случае PSQL:
psql --username ics --password ics --dbname mydb
 */
-- - 1. myschem1 для mydb1
CREATE SCHEMA IF NOT EXISTS myschem1;

-- - 2. myschem2 для mydb1
CREATE SCHEMA IF NOT EXISTS myschem2;

-- - 3. Определитб текущую схему
select current_schema();
--  current_schema
-- ----------------
--  public
-- или
SHOW search_path;
--    search_path
-- -----------------
--  "$user", public

-- - 4. Сделать myschem2 текущей
SET search_path TO myschem2;
SHOW search_path;
--  search_path
-- -------------
--  myschem2


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
Домен - опциональный констрейнт с уникальным именем в рамках схемы.
Нужны для централизованного управления атрбутами
*/
CREATE DOMAIN mydom AS
    INT CHECK (value < 7);



-- ### 9. Добавить расширение pg_freespacemap
CREATE EXTENSION pg_freespacemap;
/*
 каждый объект (например, таблица или индекс, кроме хеш-индекса) имеет Free Space Map (FSM)
 для отслеживания свободного места (таблицы).
 если filenode (файл, хранящий, наприме, таблицу) = 12345,
 то FSM хранится в файле 12345_fsm в той же директории, что и основной файл.
 FSM оргназиован как дерево FSM-страниц
 */
