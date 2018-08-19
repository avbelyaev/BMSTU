--Задание 15. Создание распределенных баз данных со связанными таблицами средствами СУБД SQL Server 2012


--1.Создать в базах данных пункта 1 задания 13 связанные таблицы.
--2.Создать необходимые элементы базы данных (представления, триггеры), 
	--обеспечивающие работу с данными связанных таблиц (выборку, вставку, изменение, удаление).


use master;
go
if DB_ID (N'lab15_1') is not null
	drop database lab15_1;
go
create database lab15_1
go


use master;
go
if DB_ID (N'lab15_2') is not null
	drop database lab15_2;
go
create database lab15_2
go


--1.Создать в базах данных пункта 1 задания 13 связанные таблицы.
--================================================
-----------------DATABASE-1-PARENT----------------
--================================================
use lab15_1;
go

if OBJECT_ID(N'writers', N'U') is not null
    drop table writers
go
create table writers 
(
    id				int				identity(1,1),
    name			varchar(35),
 
    PRIMARY KEY (id)
);
go


--================================================
-----------------DATABASE-2-CHILD-----------------
--================================================
use lab15_2;
go

if OBJECT_ID(N'books', N'U') is not null
    drop table books;
go
create table books 
(
    title           varchar(254),
    wid             int,
    rating          int,
	instock			varchar(35),
 
    PRIMARY KEY (title)
	--foreign key references не работает между 2мя бд,
	--т.ч. таблицы "связаны" только условно, через триггеры
);
go


--2.Создать необходимые элементы базы данных (представления, триггеры), 
	--обеспечивающие работу с данными связанных таблиц (выборку, вставку, изменение, удаление).
--==============================================================
--========================P-A-R-E-N-T===========================
--==============================================================

--==============================================================
-----------------------INSERT-ON-PARENT-------------------------
--==============================================================

--normal insert on parent (writers). trigger is not necessary

use lab15_1;
go
insert into writers values
	('john tolkien'),
	('george r.r. martin'),
	('andy weir');

select * from writers

--==============================================================
-----------------------DELETE-ON-PARENT-------------------------
--==============================================================

--let there be cascade delete on child table

use lab15_1;
go
if OBJECT_ID(N'writers_del', N'TR') is not null
    drop trigger writers_del
go
create trigger writers_del
	on lab15_1.dbo.writers
	instead of delete
	--after
	as
	begin
		
		delete from lab15_1.dbo.writers
			where id in (select id from deleted)

		delete from lab15_2.dbo.books
			where wid in (select id from deleted)
	end
go

delete from writers
	where writers.name like '%andy%';

select * from writers

--==============================================================
-----------------------UPDATE-ON-PARENT-------------------------
--==============================================================

--let the "id" be unupdatable 
--and the "name" be updatable so we dont care wht happens to child's 'foreign key' (wid)

use lab15_1;
go
if OBJECT_ID(N'writers_upd', N'TR') is not null
    drop trigger writers_upd
go
create trigger writers_upd
	on lab15_1.dbo.writers
	instead of update
	as
	begin

		if UPDATE(id)
			RAISERROR('[PARENT-UPD TRIGGER]: "id" cant be modified', 16, 1);

		if UPDATE(name)
		begin
			update writers
				set name = (select name from inserted)
				where writers.id = (select id from inserted)
		end

	end
go

update writers
	set name = 'george martin'
	where name like '%martin%'

select * from writers



--==============================================================
--=========================C-H-I-L-D============================
--==============================================================

--==============================================================
------------------------INSERT-ON-CHILD-------------------------
--==============================================================

--check if child's foreign key (wid) references to the existing entry in parent (writers)

use lab15_2;
go
if OBJECT_ID(N'books_ins', N'TR') is not null
    drop trigger books_ins
go
create trigger books_ins
	on lab15_2.dbo.books
	instead of insert
	as
	begin

		if exists (select wid from inserted where wid not in (select id from lab15_1.dbo.writers))
			begin
				RAISERROR('[CHILD-INS TRIGGER]: reference to parental entry was not found!', 16, 1);
			end

		else
			begin
				insert into books 
					select * from inserted
			end

	end
go

insert into books values
	('a dance with dragons', 2, 9, 'available'),
	('a dream of spring', 2, -1, 'n/a'),
	('the hobbit, or blablabla', 1, 7, 'available');

select * from books

--==============================================================
-----------------------DELETE-ON-CHILD--------------------------
--==============================================================

--just regular delete

delete from books
	where title like '%dream%';

select * from books

--==============================================================
-----------------------UPDATE-ON-CHILD--------------------------
--==============================================================

--let the "wid" and "title" be unupdatable 

use lab15_2;
go
if OBJECT_ID(N'books_upd', N'TR') is not null
    drop trigger books_upd
go
create trigger books_upd
	on lab15_2.dbo.books
	instead of update
	as
	begin

		if UPDATE(wid) or UPDATE(title)
			begin
				RAISERROR('[CHILD-UPD TRIGGER]: "wid" and "title" are unupdatable', 16, 1);
			end

		if UPDATE(rating) or UPDATE(instock)
			begin
				update books
					set
						rating = (select rating from inserted where inserted.title = books.title),
						instock = (select instock from inserted where inserted.title = books.title)
					where title = (select title from inserted where inserted.title = books.title)
			end
	end
go

update books
	set rating = 10
	where title like '%dragons%'

update books
	set instock = 'n/a'
	where title like '%hobbit%'

select * from books

