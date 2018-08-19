--Лабораторная работа №7. Представления и индекaсы

--1.Создать представление на основе одной из таблиц задания 6.
--2.Создать представление на основе полей обеих связанных таблиц задания 6.
--3.Создать индекс для одной из таблиц задания 6, включив в него дополнительные неключевые поля.
--4.Создать индексированное представление.


--run on launched lab6!
use lab6;
go

--1.Создать представление на основе одной из таблиц задания 6.
if OBJECT_ID(N'bookshelf', N'V') is not null
	drop view bookshelf;
go
create view bookshelf as 
	select *
	from books
	where year between 1930 and 1950;
go

--select * from bookshelf;
go


--2.Создать представление на основе полей обеих связанных таблиц задания 6.
if OBJECT_ID(N'bookshelf_joined', N'V') is not null
	drop view bookshelf_joined;
go
create view bookshelf_joined as
	select 
		b.name as book_name, 
		w.firstname + ' ' + w.lastname as author, 
		b.year
	from books b 
	inner join writers w
		on b.bwid = w.wid;
go

select * from bookshelf_joined;
go


--3.Создать индекс для одной из таблиц задания 6, включив в него дополнительные неключевые поля.
if EXISTS (select name from sys.indexes 
			where name = N'ix_books_bid_year')
	drop index ix_books_bid_year on books;
go
create index ix_books_bid_year
	on books (name)
	include (year);
go


--4.Создать индексированное представление.
if OBJECT_ID(N'bookshelf_ix_view', N'V') is not null
	drop view bookshelf_ix_view;
go
--table that is referenced by view_with_SCHEMABINDING cannot be modified
--if this modification affects view
create view bookshelf_ix_view 
	with SCHEMABINDING	
	as select bid, year
	from dbo.books
	where year > 1930;
go
--unique clustered index must be the !first! to be made
--only then any other index can be applied
if EXISTS (select name from sys.indexes 
			where name = N'bookshelf_ix_view')
	drop index bookshelf_ix_view on books;
go
create UNIQUE CLUSTERED index ix_books_bid_year
	on bookshelf_ix_view(bid, year);
go

