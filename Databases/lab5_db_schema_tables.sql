--Лабораторная работа №5. Операции с базой данных, файлами, схемами

--1.Создать базу данных (CREATE DATABASE…, определение настроек размеров файлов).
--2.Создать произвольную таблицу(CREATE TABLE…).
--3.Добавить файловую группу и файл данных(ALTER DATABASE…).
--4.Сделать созданную файловую группу файловой группой по умолчанию.
--5.(*) Создать еще одну произвольную таблицу.
--6.(*) Удалить созданную вручную файловую группу.
--7.Создать схему, переместить в нее одну из таблиц, удалить схему.



--1.Создать базу данных (CREATE DATABASE…, определение настроек размеров файлов).
use master;
go
if DB_ID (N'lab5') is not null
	drop database lab5;
go
create database lab5
on ( 
	NAME = lab5dat, 
	FILENAME = 'C:\Users\me\Documents\DB_Labs\lab5\lab5dat.mdf',
	SIZE = 10, 
	MAXSIZE = UNLIMITED, 
	FILEGROWTH = 5 
	)
log on ( 
	NAME = lab5log, 
	FILENAME = 'C:\Users\me\Documents\DB_Labs\lab5\lab5log.ldf',
	SIZE = 5, 
	MAXSIZE = 20, 
	FILEGROWTH = 5 
	);
go


--3.Добавить файловую группу и файл данных(ALTER DATABASE…).
use master;
go
alter database lab5
	add filegroup lab5fg;
go
alter database lab5
	add file (
		name = extrafile,
		filename = 'C:\Users\me\Documents\DB_Labs\lab5\lab5extra.ndf',
		size = 5,
		maxsize = 10,
		filegrowth = 1 
		)
	to filegroup lab5fg
go


--4.Сделать созданную файловую группу файловой группой по умолчанию.
alter database lab5
	modify filegroup lab5fg default;
go
alter database lab5
	modify filegroup [primary] default;
go


--2.Создать произвольную таблицу(CREATE TABLE…).
use lab5;
go
if OBJECT_ID(N'owner', N'U') is not null
	drop table owner;
go
create table owner (
	oid 		int not 	null,
	ownername 	varchar(35) 	not null,
	groupname 	varchar(35) 	not null,
	privelegelvl 	int 		not null,
	email 		varchar(254) 	null,
	ostatus 	bit 		not null
	);
go

insert into owner(oid, ownername, groupname, privelegelvl, email, ostatus)
	values (0, 'me', 'Users', 100500, 'me@example.com', 1);

select * from owner;
go


--5.(*) Создать еще одну произвольную таблицу.
if OBJECT_ID(N'extratable', N'U') is not null
	drop table extratable;
go
create table extratable (
	name 	varchar(35) 	null
	);
go


--6.(*) Удалить созданную вручную файловую группу.
--сначала удалить из нее файл
alter database lab5
	remove file extrafile;
go
alter database lab5
	remove filegroup lab5fg;
go



--7.Создать схему, переместить в нее одну из таблиц, удалить схему.
if SCHEMA_ID(N'lab5schema') is not null
	drop schema lab5schema;
go
create schema lab5schema;
go

--modify schema
alter schema lab5schema 
	transfer extratable;
go
if OBJECT_ID(N'lab5schema.extratable', N'U') is not null
	drop table lab5schema.extratable;
go

--drop schema
drop schema lab5schema;
go

