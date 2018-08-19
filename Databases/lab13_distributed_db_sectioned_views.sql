--Задание 13. Создание распределенных баз данных на основе секционированных представлений


--1.Создать две базы данных на одном экземпляре СУБД SQL Server 2012.
--2.Создать в базах данных п.1. горизонтально фрагментированные таблицы.
--3.Создать секционированные представления, обеспечивающие работу с данными таблиц (выборку, вставку, изменение, удаление).



--1.Создать две базы данных на одном экземпляре СУБД SQL Server 2012.
use master;
go
if DB_ID (N'lab13_1') is not null
	drop database lab13_1;
go
create database lab13_1
go


use master;
go
if DB_ID (N'lab13_2') is not null
	drop database lab13_2;
go
create database lab13_2
go



--2.Создать в базах данных п.1. горизонтально фрагментированные таблицы.
use lab13_1;
go
if OBJECT_ID(N'dbo.books', N'U') is not null
    drop table dbo.books;
go
create table dbo.books (
	id				int not null,
    title           varchar(254),
    writer			varchar(35),
 
    PRIMARY KEY (id),
	CONSTRAINT CHK_books_id
                CHECK (id <= 5)
    );
go


use lab13_2;
go
if OBJECT_ID(N'dbo.books', N'U') is not null
    drop table dbo.books;
go
create table dbo.books (
	id				int not null,
    title           varchar(254),
    writer			varchar(35),
 
    PRIMARY KEY (id),
	CONSTRAINT CHK_books_id
                CHECK (id > 5)
    );
go



--3.Создать секционированные представления, обеспечивающие работу с данными таблиц 
	--(выборку, вставку, изменение, удаление).
use lab13_1;
go
if OBJECT_ID(N'horizontal_dist_v', N'V') is not null
	drop view horizontal_dist_v;
go
create view horizontal_dist_v as
	select * from lab13_1.dbo.books
	union all					
	select * from lab13_2.dbo.books
go



insert horizontal_dist_v values
	(1, 'the lord of the rings', 'john tolkien'),
	(2, 'the hobbit, or there and back again', 'john tolkien'),
	(3, 'the winds of winter', 'unknown'),
	(6, 'a game of thrones', 'george martin'),
	(10, 'a dream of spring', 'John Tolkien');

select * from horizontal_dist_v

update horizontal_dist_v
	set writer = 'george martin'
	where title like '%dream%'

delete horizontal_dist_v
	where writer = 'unknown'

select * from lab13_1.dbo.books;
select * from lab13_2.dbo.books;

