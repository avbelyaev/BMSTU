
### 1. Сервер

#### 1. Запустить сервер.
`docker-compose up`

#### 2. Остановить сервер
Ctrl-C

#### 3. Посмотреть журнал.
журнал - в `./postgres-logs`

пример записей:
```
2019-10-08 23:07:19.192 UTC [24] LOG:  database system was not properly shut down; automatic recovery in progress
2019-10-08 23:07:19.289 UTC [24] LOG:  invalid record length at 0/1638080: wanted 24, got 0
2019-10-08 23:07:19.289 UTC [24] LOG:  redo is not required
2019-10-08 23:07:19.508 UTC [1] LOG:  database system is ready to accept connections
2019-10-08 23:16:46.606 UTC [61] ERROR:  database "ics" already exists
2019-10-08 23:16:46.606 UTC [61] STATEMENT:  CREATE DATABASE ics
        WITH 
        OWNER = ics
        ENCODING = 'UTF8'
        LC_COLLATE = 'en_US.utf8'
        LC_CTYPE = 'en_US.utf8'
        TABLESPACE = myts1
        CONNECTION LIMIT = -1;
2019-10-08 23:17:26.265 UTC [61] ERROR:  cannot drop the currently open database
2019-10-08 23:17:26.265 UTC [61] STATEMENT:  DROP DATABASE ics;
```

#### 4. Определить версию сервера
```postgres-sql
SELECT verison();
-- PostgreSQL 10.10 (Debian 10.10-1.pgdg90+1) 
--   on x86_64-pc-linux-gnu, compiled by gcc
```




### 2. Настройки сервера

#### 1. Посмотрите файл настройки сервера
в контейнере под адресу `/var/lib/postgresql/data/postgresql.conf`

- аутнетификация и авторзация
- файлы
- WAL
- соединение
- оптимизация 
- обработка ошибок
- vacuum

#### 2. Посмотреть содержание pg_log, pg_сlog, pg_хlog

в контейнере под адресу `/var/lib/postgresql/data`

- pg_log - является по умолчанию местом, где хранятся журналы деятельности. 
Они включают в себя сообщения об ошибках, записи о запросах, 
и сообщения во время старта\выключения СУБД.

- pg_xlog (PG10 -> pg_wal) — это место, где PostgreSQL хранит журнал транзакций. 
Этот набор бинарных файлов, с названиями вида '00000001000000000000008E', 
которые содержат образы данных последних транзакций. 

- pg_clog (PG10 -> pg_xact) - содержит журналы метаданных транзакций. 
Этот журнал говорит серверу, какие транзакции завершены, а какие нет. 
Если вы когда-нибудь удалите файлы из pg_clog, вы можете смело удалить и весь каталог 
базы данных. Не существует способа восстановить базу данных без этих журналов.

#### 3. Посмотреть настройки конфигурации 
посмотрели

#### 4. Посмотреть значение параметра shared_buffers

`shared_buffers = 128MB			# min 128kB`

#### 5. Увеличите в три раза 'heap memory per user' процедурой sp_configure

какой нафиг sp_configure. это вообще SQL server

#### 6. Изменить параметр файла конфигурации

Изменили `shared_buffers = 256MB` в `custom-pg.conf`. Рестартанули контейнер




### 3. Создать tablespace

#### 1. Создать tablespace “myts1” (myts2) в директории MyDB1 (MyDB2) с пользователем postgres

- в контейнере:
 
```bash
mkdir /home/space1 && chown postgres /home/space1
mkdir /home/space2 && chown postgres /home/space2
```

- по отдельности запускаем:

```postgres-sql
CREATE TABLESPACE myts1 OWNER ics LOCATION '/home/space1';
CREATE TABLESPACE myts2 OWNER ics LOCATION '/home/space2';
```



### 4. Создать БД

#### 1. Создать БД mydb

```postgres-sql
DROP DATABASE IF EXISTS mydb;

CREATE DATABASE mydb
    WITH 
    OWNER = ics
    ENCODING = 'UTF8'
    TABLESPACE = myts1;
```

#### 2. Создать БД mytest

```postgres-sql
CREATE DATABASE mytest
    WITH 
    OWNER = ics
    ENCODING = 'UTF8'
    TABLESPACE = myts2;
```

#### 3. Уничтожить БД mytest

```postgres-sql
DROP DATABASE IF EXISTS mytest;
```



### 5. Создать схемы

#### 1. myschem1 (myschem2) для mydb1

- переключаемся на БД mydb1
```postgres-sql
CREATE SCHEMA IF NOT EXISTS myschem1;
CREATE SCHEMA IF NOT EXISTS myschem2;
```

#### 3. Определит текущею схему

`select current_schema();` или `SHOW search_path;`

#### 4. Сделать myschem2 текущей

```postgres-sql
SET search_path TO myschem2;
SHOW search_path;
```



### 6. Создать последовательность

#### 1. MySq c шагом 1

```postgres-sql
CREATE SEQUENCE IF NOT EXISTS mysq1
    INCREMENT 1
    START 1;
```


### 7. Создать Тип

#### 1.Создать новый тип с описанием работника

The CREATE TYPE statement allows you to create a composite type, 
which can be use as the return type of a function.

```postgres-sql
CREATE TYPE fio AS
(
    name    character(40), 
    soname  character(40), 
    family  character(40),
    gender  character(1)
);
```



### 8. Создать домен

In PostgreSQL, a domain is a data type with optional constraints e.g., 
NOT NULL, CHECK etc. A domain has a unique name within the schema scope.
Domains are useful for centralizing management of fields with the common constraints.

```postgres-sql
CREATE DOMAIN mydom AS 
    INT CHECK (value < 7);
```



### 9. Добавить расширение pg_freespacemap

`CREATE EXTENSION pg_freespacemap;`



### 10.Таблицы

- a. Создать в схеме myschem1 таблицы согласно перечисленным ниже требованиям
- b. Создать связи между таблицами, как показано на логической схеме БД
- c. Посмотреть каталог где располагается таблица


- Таблица EMPLOYERS
    - 1. Первичные ключи объявить на основе последовательности.
    - 2. aSTelephone объявить как массив строк.
    - 3. Добавить столбец FIO1 типа fio

- Таблица требований TIT_OUT
    - 1. Первичные ключи объявить на основе автоинкрементальный (Serial).
    - 2. iCODE уникальное поле

- Таблица NOM_OUT
    1. Первичные ключи определиться как автоинкрементальный (Serial)
    2. Дата по умолчанию текущая
    3. Предусмотреть поля для внешних ключей ID_TIT и ID_NOM
    4. Связать таблицы как указано на схеме

- Таблица TIT_IN
    1. Первичные ключи объявить автоинкрементальный
    2. iCODE уникальный
    
- Таблицы NOM_IN
    1. Первичные ключи определиться как автоинкрементальный (Serial).
    2. Предусмотреть поля для внешних ключей ID_TIT и ID_NOM
    3. Дата текущая
    4. Связать таблицы как указано на схеме

- Таблица NOMECLATURA
    1. Первичные ключи объявить как автоинкрементальный (Serial)..
    2. SCODE unique
    3. dtINPOUT текущая дата  // такого вообще нет
    
- Таблица HEAD
    1. Наследованием EMPLOYES 
    2. Дополнительное поле nadd - надбавка
    
- Таблица INSTRUCTION
    1. поле для текстовой информации inst
    2. поле ins поле типа tsvector
    
the tsvector data type, where ts stands for "text search"); 
to_tsquery for querying the vector for occurrences of certain words or phrases.


    
```postgres-sql
DROP TABLE IF EXISTS head;
DROP TABLE IF EXISTS employers;
DROP TABLE IF EXISTS nom_in;
DROP TABLE IF EXISTS nom_out;
DROP TABLE IF EXISTS tit_out;
DROP TABLE IF EXISTS tit_in;
DROP TABLE IF EXISTS nomenclatura;
DROP TABLE IF EXISTS instruction;



CREATE TABLE employers
(
    ID_EMP              BIGINT PRIMARY KEY DEFAULT nextval('mysq1'),
    SName               VARCHAR(40),
    SFamily             VARCHAR(40),
    FIO1                fio,
    SPosition           VARCHAR(40),
    SSex                BOOLEAN,
    aSTelephone         TEXT [],
    SOrganization       VARCHAR(255),
    S_FIO_EMPL          VARCHAR(50),
    SLogin              VARCHAR(40),
    SPSW                CHAR(10),
    IAccess             BOOLEAN,
    sTCHF               VARCHAR(20)
);


CREATE TABLE nomenclatura
(
    ID_DRG               SERIAL PRIMARY KEY,
    REM_ID               NUMERIC(6),
    sCODE                VARCHAR(20) UNIQUE,
    DRUG_BARF            NUMERIC(14),
    DRUG_NAME            VARCHAR(255),
    PACK1_QTTY           NUMERIC(4),
    NOM_QTTY             NUMERIC(4),
    INVAL_DATE           DATE,
    siACTUAL             SMALLINT,
    sINSTRUCTION         TEXT,
    CHECK_DATE           DATE,
    sNOTE                TEXT,
    sTCHF                VARCHAR(20)
);

CREATE TABLE tit_out
(
    ID                   SERIAL PRIMARY KEY,
    iDIVISION            INTEGER,
    iCODE                NUMERIC(6) UNIQUE,
    dDATE                DATE,
    dLAST_CORR           DATE,
    iWORK_UP             INT,
    dWORK_UP             TIMESTAMP WITHOUT TIME ZONE,
    sSUBSCR_APT          VARCHAR(30),
    sSUBSCR_DIV          VARCHAR(30),
    sNOTE                TEXT,
    sTCHF                VARCHAR(20)
);

CREATE TABLE tit_in
(
    ID_TIT               SERIAL PRIMARY KEY,
    EMP_ID               NUMERIC(10),
    iCODE                VARCHAR(16) UNIQUE,
    dDATE                DATE,
    dLASTCORR            DATE,
    dWORKUP              TIMESTAMP WITHOUT TIME ZONE,
    nSUM                 NUMERIC(11, 2),
    nCORR                NUMERIC(11, 2),
    sSUBSCR              VARCHAR(30),
    sNOTE                TEXT,
    sTCHF                VARCHAR(20)
);


CREATE TABLE nom_in
(
    ID_CL                SERIAL PRIMARY KEY,
    nQUANTITY            NUMERIC(12, 3),
    nQUANT               NUMERIC(12, 3),
    dLASTCORR            DATE NOT NULL DEFAULT CURRENT_DATE,
    dtINPUT              TIMESTAMP WITHOUT TIME ZONE,
    dDATE                DATE NOT NULL DEFAULT CURRENT_DATE,
    PACK1_QTTY           NUMERIC(4),
    NOM_QTTY             NUMERIC(4),
    siERROR              SMALLINT,
    sTCHF                VARCHAR(20),
    ID_TIT               INT,   --nullable
    ID_NOM               INT,   --nullable
    CONSTRAINT fk_tit_out 
        FOREIGN KEY (ID_TIT) 
        REFERENCES tit_out (ID),
    CONSTRAINT fk_nomenclatura
        FOREIGN KEY (ID_NOM)
        REFERENCES nomenclatura (ID_DRG)
);

CREATE TABLE nom_out
(
    ID                   SERIAL PRIMARY KEY,
    nQUANT_C             NUMERIC(12, 3),
    nQUANT_F             NUMERIC(12, 3),
    dLAST_COOR           DATE NOT NULL DEFAULT CURRENT_DATE,
    siERROR              SMALLINT,
    dtINPUT              TIMESTAMP WITHOUT TIME ZONE,
    dDATE                DATE NOT NULL DEFAULT CURRENT_DATE,
    sUNIT                VARCHAR(1),
    nPRICE               NUMERIC(12, 2),
    siPROPIS             SMALLINT,
    TIT_ID               NUMERIC(10),
    sTCHF                VARCHAR(20),
    ID_TIT               INT,    --nullable
    ID_NOM               INT,    --nullable
    CONSTRAINT fk_tit_out 
        FOREIGN KEY (ID_TIT) 
        REFERENCES tit_out (ID),
    CONSTRAINT fk_nomenclatura
        FOREIGN KEY (ID_NOM)
        REFERENCES nomenclatura (ID_DRG)
);

CREATE TABLE head
(
	nadd                INT,
	PRIMARY KEY (ID_EMP)
) 
INHERITS (employers);

CREATE TABLE instruction
(
	inst                TEXT,
	ins                 tsvector
);
```


- Временная таблица Temp EMPLOYERS
    1. Структура таблицы такая же как и EMPLOYERS
    
Временная таблица существует только в рамках текущей сессии. 
Если создать новое подключение и выплнить `SELECT * FROM temp_employers;`, то таблица не найдется
    
```postgres-sql
DROP TABLE IF EXISTS temp_employers;

CREATE TEMPORARY TABLE temp_employers
(
    ID_EMP              BIGINT PRIMARY KEY DEFAULT nextval('mysq1'),
    SName               VARCHAR(40),
    SFamily             VARCHAR(40),
    FIO1                fio,
    SPosition           VARCHAR(40),
    SSex                BOOLEAN,
    aSTelephone         TEXT [],
    SOrganization       VARCHAR(255),
    S_FIO_EMPL          VARCHAR(50),
    SLogin              VARCHAR(40),
    SPSW                CHAR(10),
    IAccess             BOOLEAN,
    sTCHF               VARCHAR(20)
);

SELECT * FROM temp_employers;
```


Создать индекс по коду накладной, и дате поступления.

```postgres-sql
CREATE INDEX idx_nomanclatura_scode_inval_date
    ON nomenclatura(sCODE, INVAL_DATE);
```
