--Лабораторная работа №6. Ключи, ограничения, значения по умолчанию

--1.Создать таблицу с автоинкрементным первичным ключом.
--2.Добавить поля, для которых используются ограничения (CHECK), значения по умолчанию(DEFAULT), также использовать функции для вычисления.
--3.Создать таблицу с первичным ключом на основе глобального уникального идентификатора.
--4.Создать таблицу с первичным ключом на основе последовательности.
--5.Создать две связанные таблицы, и протестировать на них различные варианты действий для ограничений ссылочной целостности (NO ACTION| CASCADE | SET NULL | SET DEFAULT).

--create db
use master;
go
if DB_ID (N'lab6') is null
	create database lab6
	on ( 
		NAME = lab6dat, 
		FILENAME = 'C:\Users\me\Documents\DB_Labs\lab6\lab6dat.mdf',
		SIZE = 10, 
		MAXSIZE = UNLIMITED, 
		FILEGROWTH = 5 
		)
	log on ( 
		NAME = lab6log, 
		FILENAME = 'C:\Users\me\Documents\DB_Labs\lab6\lab6log.ldf',
		SIZE = 5, 
		MAXSIZE = 20, 
		FILEGROWTH = 5 
		);
go


use lab6;
go
--1.Создать таблицу с автоинкрементным первичным ключом.
--2.Добавить поля, для которых используются ограничения (CHECK), значения по умолчанию(DEFAULT), также использовать функции для вычисления.
if OBJECT_ID(N'owners', N'U') is not null
	drop table owners
go
create table owners (
	oid			int IDENTITY(100,100)		not null,
	firstname		varchar(35)			not null,
	lastname		varchar(35)			null,
	ostatus			int				null
		CONSTRAINT DF_owners_status DEFAULT (1),
	controlvalue		int				not null,
	createdate		datetime			not null
		CONSTRAINT DF_owners_createdate DEFAULT (getdate()),

	PRIMARY KEY (oid),
	CONSTRAINT CHK_owners_controlvalue
		CHECK (controlvalue > 500)
	);
go

insert into owners(firstname, controlvalue)
	values 
	('me', 3405),
	('Andrew', 1203), 
	('Johan', 501);
	
go
select * from owners;

select IDENT_CURRENT ('owners') AS current_id;
--300

insert into owners(firstname, controlvalue) 
	values
		('Rico', 5553535),
		('Sam', 8800),
		('Jo', 3222),
		('Q', 1007);
go


--3.Создать таблицу с первичным ключом на основе глобального уникального идентификатора.
if OBJECT_ID(N'usergroups', N'U') is not null
	drop table usergroups
go
create table usergroups (
	gid		UNIQUEIDENTIFIER DEFAULT NEWID(),
	groupname	varchar(35),
	usersnumber 	int
		CONSTRAINT DF_usergroups_usersnumber DEFAULT (0),

	PRIMARY KEY (gid),
	CONSTRAINT CHK_usergroups_usersnumber
		CHECK (usersnumber >= 0)
	);
go

insert into usergroups(groupname, usersnumber)
	values 
		('admins', 1),
		('stdandardusers', 100500);
go
insert into usergroups values (NEWID(), 'superusers', 2);

select * from usergroups;


--4.Создать таблицу с первичным ключом на основе последовательности.
if OBJECT_ID(N'math', N'U') is not null
	drop table math
go
create table math (
	multiplication	varchar(5)	not null,
	firstnum	varchar(5)	not null,
	secondnum	int		not null,

	PRIMARY KEY (multiplication)
	);
go


--drop SEQUENCE lab6schema.incmul;
if EXISTS(select * from sys.objects where object_id = OBJECT_ID(N'lab6schema.incmul') and type = 'SO')
	drop SEQUENCE lab6schema.incmul;
		
if SCHEMA_ID(N'lab6schema') is not null
	drop schema lab6schema;
go
create schema lab6schema;
go

create SEQUENCE lab6schema.incmul
	START WITH 0
	INCREMENT BY 2
	MAXVALUE 10;
go

insert math values 
		('0 = ', '0 * ', NEXT VALUE FOR lab6schema.incmul),
		('2 = ', '1 * ', NEXT VALUE FOR lab6schema.incmul),
		('8 = ', '2 * ', NEXT VALUE FOR lab6schema.incmul);

select * from math;





use lab6;
go 
--5.Создать две связанные таблицы, и протестировать на них различные варианты действий 
--для ограничений ссылочной целостности (NO ACTION| CASCADE | SET NULL | SET DEFAULT).
if (OBJECT_ID(N'FK_books_bwid', N'F') is not null)
begin
	alter table books
		drop CONSTRAINT FK_books_bwid
end
go
/*C = CHECK constraint
D = DEFAULT (constraint or stand-alone)
F = FOREIGN KEY constraint
PK = PRIMARY KEY constraint
UQ = UNIQUE constraint*/

if OBJECT_ID(N'writers', N'U') is not null
	drop table writers
go
create table writers (
	wid		int			not null,
	firstname	varchar(35)		not null,
	midname		varchar(35)		null,
	lastname	varchar(35)		not null,
	tmp		int

	PRIMARY KEY (wid)
	);
go

insert writers values 
	(1, 'john', 'r.r.', 'tolkien', 1),
	(2, 'george', 'r.r.', 'martin', 2),
	(3, 'herbert', 'g.', 'wells', 3),
	(4, 'george', '', 'orwell', 4),
	(5, 'raymond', 'd.', 'bradbury', 5);
go

select * from writers;


--=====from lab7
--if OBJECT_ID(N'bookshelf_ix_view', N'V') is not null
--	drop view bookshelf_ix_view;
--go
--=====
if OBJECT_ID(N'books', N'U') is not null
	drop table books
go
create table books (
	bid	int IDENTITY(1, 1)	not null,
	name	varchar(254)		not null,
	bwid	int			not null
		CONSTRAINT DF_books_bwid DEFAULT (3),
	year	int			null,

	PRIMARY KEY (bid),
	CONSTRAINT FK_books_bwid
		FOREIGN KEY (bwid) 
		REFERENCES writers(wid) 
		ON DELETE SET DEFAULT
	);
go

insert into books(name, bwid, year)
	values
		('the lord of the rings', 1, 1954),
		('the hobbit, or there and back again', 1, 1937),
		('a song of ice and fire', 2, 1996),
		('1984', 4, 1949),
		('fahrenheit 451', 5, 1953);
go

select * from books;


/*
	NO ACTION -> Error
	SET DEAFULT -> set VALID deafult
	SET NULL -> -||- (CARE OF "NOT_NULL" value in books.bwid)
	CASCADE (del) -> -||-
*/

--delete from writers where writers.tmp = 1;
--go

select * from writers;
select * from books;

