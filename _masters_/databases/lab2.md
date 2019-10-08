
### 1. Сервер

#### 1. Запустить сервер.
`docker-compose up`

#### 2. Остановить сервер
Ctrl-C

#### 3. Посмотреть журнал.
журнал - в `./postgres-logs`

#### 4. Определить версию сервера
```postgres-sql
SELECT verison();
-- PostgreSQL 10.10 (Debian 10.10-1.pgdg90+1) 
--   on x86_64-pc-linux-gnu, compiled by gcc
```




### 2. Настройки сервера

#### 1. Посмотрите файл настройки сервера
в контейнере под адресу `/var/lib/postgresql/data/postgresql.conf`

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
	name character(40), 
 	soname character(40), 
	family character(40),
	gender character(1)
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


#### Таблица EMPLOYERS

- 1. Первичные ключи объявить на основе последовательности.
- 2. aSTelephone объявить как массив строк.
- 3. Добавить столбец FIO1 типа fio

```postgres-sql
CREATE TABLE employers
(
	ID_EMP 				BIGINT PRIMARY KEY DEFAULT nextval('mysq1'),
	SName 				VARCHAR(40),
	SFamily 			VARCHAR(40),
	SPosition 			VARCHAR(40),
	SSex 				BOOLEAN,
	aSTelephone 		TEXT [],
	SOrganization 		VARCHAR(255),
	S_FIO_EMPL 			VARCHAR(50),
	SLogin				VARCHAR(40),
	SPSW 				CHAR(10),
	IAccess 			BOOLEAN,
	sTCHF				VARCHAR(20)
)
```


#### Таблицы требований TIT_OUT

- 1. Первичные ключи объявить на основе автоинкрементальный (Serial).
- 2. iCODE уникальное поле
