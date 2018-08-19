--Лабораторная работа №9. Триггеры DML.


--1.Для одной из таблиц пункта 2 задания 7 создать триггеры на вставку, удаление и добавление, 
	--при выполнении заданных условий один из триггеров должен инициировать возникновение ошибки (RAISERROR / THROW).
--2.Для представления пункта 2 задания 7 создать триггеры на вставку, удаление и добавление, 
	--обеспечивающие возможность выполнения операций с данными непосредственно через представление.



use master;
go
if DB_ID (N'lab9') is null
	create database lab9
		on (
			NAME = lab9dat,
			FILENAME = 'C:\Users\me\Documents\DB_Labs\lab9\lab9dat.mdf',
			SIZE = 5, MAXSIZE = UNLIMITED, FILEGROWTH = 5 )
        log on (
			NAME = lab8log,
			FILENAME = 'C:\Users\me\Documents\DB_Labs\lab9\lab9log.ldf',
			SIZE = 5, MAXSIZE = 20, FILEGROWTH = 5 );
go
 
use lab9;
go


if (OBJECT_ID(N'FK_books_wid', N'F') is not null)
	alter table dbo.books
		drop CONSTRAINT FK_books_wid
go
if OBJECT_ID(N'dbo.writers', N'U') is not null
	drop table dbo.writers
go
--parental table
create table dbo.writers (
	id				int identity(1,1),
	name			varchar(35),

	PRIMARY KEY (id)
	);
go


if OBJECT_ID(N'dbo.books', N'U') is not null
	drop table dbo.books;
go
--child table
create table dbo.books (
	title			varchar(254),
	wid				int
		CONSTRAINT DF_books_wid DEFAULT(0),
	year			int,
	rating			int
		CONSTRAINT DF_books_rating DEFAULT(0),

	PRIMARY KEY (title),
	CONSTRAINT FK_books_wid
		FOREIGN KEY (wid)
		REFERENCES dbo.writers(id)
	);
go


if OBJECT_ID(N'dbo.joined_view', N'V') is not null
	drop view dbo.joined_view;
go
create view dbo.joined_view as
	select 
		b.title			as title,
		w.name			as author, 
		b.year			as year,
		b.rating		as rating
	from dbo.books b 
	inner join dbo.writers w
		on b.wid = w.id;
go

insert dbo.writers values
	('john tolkien'),
	('george martin'),
	('george orwell'),
	('andy weir'),
	('raymond bradbury'),
	('herbert wells'),
	('daniel brown');
go

insert dbo.books values
	('the lord of the rings', 1, 1954, 10),
	('the hobbit, or there and back again', 1, 1937, 8),
	('a game of thrones', 2, 1996, 7),
	('a dance with dragons', 2, 2011, 9),
	('1984', 3, 1949, 6),
	('the martian', 4, 2011, 8),
	('fahrenheit 451', 5, 1953, 6),
	('the war of the worlds', 6, 1897, 7),
	('the da vinci code', 7, 2003, 6);
go
--this query ^ is used to get tables back to their normal states





//ANOTHER QUERY:



use lab9;
go

if OBJECT_ID(N'dbo.rand_insert', N'P') is not null
	drop procedure dbo.rand_insert
go
CREATE PROCEDURE dbo.rand_insert
AS
	SET NOCOUNT ON;
	declare @rand_ptr int = ABS(Checksum(NewID()) % 3);

	if (@rand_ptr = 0)
		INSERT INTO dbo.writers values ('ivan ivanov');
	if (@rand_ptr = 1)
		INSERT INTO dbo.writers values ('petr petrov');
	if (@rand_ptr = 2)
		INSERT INTO dbo.writers values ('sergey sergeev');
go


--==============================================================
------------------------SIMPLE-TRIGGERS-------------------------
--==============================================================
--1.Для одной из таблиц пункта 2 задания 7 создать триггеры на вставку, удаление и добавление, 
	--при выполнении заданных условий один из триггеров должен инициировать возникновение ошибки (RAISERROR / THROW).

--trigger on INSERT with RAISERROR 
if OBJECT_ID(N'dbo.writers_trigger_insert', N'TR') is not null
	drop trigger dbo.writers_trigger_insert
go
create trigger dbo.writers_trigger_insert
	on dbo.writers
	for insert
	as
		declare @min_val int = 10;

		if exists (select *
			from inserted
			where inserted.id >= @min_val)
			RAISERROR('[INS/UPD TRIGGER]: Entry with "writer.id" > 10 were added', 10, 1);
go


--trigger on UPDATE
if OBJECT_ID(N'dbo.writers_trigger_update', N'TR') is not null
	drop trigger dbo.writers_trigger_update
go
create trigger dbo.writers_trigger_update
	on dbo.writers
	for update
	as
		print 'table dbo.writers has been updated'
go


--trigger on DELETE
if OBJECT_ID(N'dbo.writers_trigger_delete', N'TR') is not null
	drop trigger dbo.writers_trigger_delete
go
create trigger dbo.writers_trigger_delete
	on dbo.writers
	instead of delete
	as
		declare @comparator int = 8;
		--delete only if there are more than 10 entries
		if (select count(*) from dbo.writers) > 10
		begin
			print '[DEL TRIGGER]: Entries with "w_id" >= ' + CAST(@comparator as varchar) + ' are deleted'
			delete from dbo.writers
				where writers.id >= @comparator
		end
go
--EXEC dbo.rand_insert;
--delete from dbo.writers where id > 1;
--select * from dbo.writers;



--==============================================================
-------------------------INSERT-ON-VIEW-------------------------
--==============================================================
--2.Для представления пункта 2 задания 7 создать триггеры на вставку, удаление и добавление, 
	--обеспечивающие возможность выполнения операций с данными непосредственно через представление.

--simple triggers modified base table so
	--rerun first query (create tables and view)
	--and disable triggers above
disable trigger dbo.writers_trigger_delete on dbo.writers;
go
disable trigger dbo.writers_trigger_update on dbo.writers;
go
disable trigger dbo.writers_trigger_insert on dbo.writers;
go

if OBJECT_ID(N'dbo.joined_trigger_insert', N'TR') is not null
	drop trigger dbo.joined_trigger_insert
go
create trigger dbo.joined_trigger_insert
	on dbo.joined_view
	instead of insert
	as
	begin

		insert into dbo.writers
			select distinct i.author
				from inserted as i
				where i.author not in (select name
					from dbo.writers)

		insert into dbo.books
			select 
					i.title, 
					(select id from dbo.writers as w where i.author = w.name), 
					i.year, 
					i.rating
				from inserted as i
	end
go


insert into dbo.joined_view values
	('the invisible man', 'herbert wells', 1897, 6)
insert into dbo.joined_view values
	('a dream of spring', 'george martin', 2016, 9)
insert into dbo.joined_view values
	('RANDOM BOOK', 'NEW AUTHOR', 111, 10)
select * from dbo.joined_view
select * from dbo.books
select * from dbo.writers



--==============================================================
-------------------------DELETE-ON-VIEW-------------------------
--==============================================================
if OBJECT_ID(N'dbo.joined_trigger_delete', N'TR') is not null
	drop trigger dbo.joined_trigger_delete
go
create trigger dbo.joined_trigger_delete
	on dbo.joined_view
	instead of delete
	as
	begin
		delete from dbo.books
			where books.title in (select d.title
				from deleted as d)
	end
go


delete from dbo.joined_view
	where joined_view.author in ('daniel brown', 'andy weir')
delete from dbo.joined_view
	where joined_view.title = 'the invisible man'
delete from dbo.joined_view
	where joined_view.rating = 10
select * from dbo.joined_view
select * from dbo.books



--==============================================================
-------------------------UPDATE-ON-VIEW-------------------------
--==============================================================
if OBJECT_ID(N'dbo.joined_trigger_update', N'TR') is not null
	drop trigger dbo.joined_trigger_update
go
create trigger dbo.joined_trigger_update
	on dbo.joined_view
	instead of update
	as 
	begin

		if UPDATE(title) or UPDATE(author) 
			RAISERROR('[UPD TRIGGER]: "title" and "author" cant be modified', 16, 1)

		if UPDATE(year) or UPDATE(rating)
			update dbo.books
				set 
					books.year = (select year from inserted where inserted.title = books.title),
					books.rating = (select rating from inserted where inserted.title = books.title)
				where books.title = (select title from inserted where inserted.title = books.title)

	end
go

update dbo.joined_view 
	set rating = 451, 
		year = 451
	where joined_view.author = 'raymond bradbury';
update dbo.joined_view
	set rating = 10
	where joined_view.author = 'george martin';
update dbo.joined_view 
	set year = 2045 
	where joined_view.title = 'a dream of spring';
select * from dbo.joined_view
select * from dbo.books

