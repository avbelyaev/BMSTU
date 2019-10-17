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
