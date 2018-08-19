--Лабораторная работа №8. Хранимые процедуры, курсоры и пользовательские функции


--1.Создать хранимую процедуру, производящую выборку из некоторой таблицы и возвращающую результат выборки в виде курсора.
--2.Модифицировать хранимую процедуру п.1. таким образом, чтобы выборка осуществлялась с формированием столбца, значение которого формируется пользовательской функцией.
--3.Создать хранимую процедуру, вызывающую процедуру п.1., осуществляющую прокрутку возвращаемого курсора и выводящую сообщения, 
    --сформированные из записей при выполнении условия, заданного еще одной пользовательской функцией.
--4.Модифицировать хранимую процедуру п.2. таким образом, чтобы выборка формировалась с помощью табличной функции.


use master;
go
if DB_ID (N'lab8') is null
        create database lab8
        on (
                NAME = lab8dat,
                FILENAME = 'C:\Users\me\Documents\DB_Labs\lab8\lab8dat.mdf',
                SIZE = 10,
                MAXSIZE = UNLIMITED,
                FILEGROWTH = 5
                )
        log on (
                NAME = lab8log,
                FILENAME = 'C:\Users\me\Documents\DB_Labs\lab8\lab8log.ldf',
                SIZE = 5,
                MAXSIZE = 20,
                FILEGROWTH = 5
                );
go
 
use lab8;
go
if OBJECT_ID(N'dbo.books', N'U') is not null
	drop table dbo.books;
go
create table dbo.books (
	bid			int			not null	IDENTITY(1, 1),
	book_name		varchar(254)		not null,
	author_lastname		varchar(35)		not null,
	author_initials		varchar(35)		null,
	year			int			not null,
	genre			varchar(254)		null,

	PRIMARY KEY (bid)
);
go

insert into dbo.books(book_name, author_lastname, author_initials, genre, year)
	values
		('the lord of the rings', 'tolkien', 'john r. r.', 'high fantasy', 1954),
		('the hobbit, or there and back again', 'tolkien', 'john r. r.', 'high fantasy', 1937),
		('a song of ice and fire', 'martin', 'george r. r.', 'epic fantasy', 1996),
		('1984', 'orwell', 'george', 'dystopian', 1949),
		('the martian', 'weir', 'andy', 'science fiction', 2011),
		('fahrenheit 451', 'bradbury', 'raymond d.', 'dystopian', 1953),
		('the time machine', 'wells', 'herbert g.', 'science fiction', 1895),
		('the war of the worlds', 'wells', 'herbert g.', 'science fiction', 1897),
		('the da vinci code', 'brown', 'daniel', 'mystery-detective', 2003);
go

select * from books;
go


--==============================================================
--------------------------FUNCTIONS-----------------------------
--==============================================================
if OBJECT_ID(N'dbo.get_antiquity', N'FN') is not null
	drop function dbo.get_antiquity;
go
create function dbo.get_antiquity(@release_year int)
	returns int
	with execute as caller
	as
	begin
		declare @current_date datetime = GETDATE();
		declare @current_year int, @year_delta int;

		set @current_year = YEAR(@current_date);
		set @year_delta = @current_year - @release_year;
	
		return @year_delta;
	end
go


if OBJECT_ID(N'dbo.compare_years', N'FN') is not null
	drop function dbo.compare_years;
go
create function dbo.compare_years(@compare_what int, @compare_to int)
	returns int
	with execute as caller
	as
	begin
		declare @retval int;

		if (@compare_what >= @compare_to)
			set @retval = 1;
		else 
			set @retval = 0;

		return @retval;
	end
go


--==============================================================
--------------------------PROCEDURES----------------------------
--==============================================================
--Создать хранимую процедуру, производящую выборку из некоторой таблицы и возвращающую результат выборки в виде курсора.

--Модифицировать хранимую процедуру п.1. таким образом, чтобы выборка осуществлялась 
--с формированием столбца, значение которого формируется пользовательской функцией.

if OBJECT_ID(N'dbo.sub_proc', N'P') is not null
	drop procedure dbo.sub_proc
go
CREATE PROCEDURE dbo.sub_proc
	@curs CURSOR VARYING OUTPUT
AS
	SET NOCOUNT ON;
	SET @curs = CURSOR
	SCROLL STATIC FOR				
		select book_name, author_lastname, dbo.get_antiquity(year)  --custom column
		from dbo.books
		OPTION (MAXRECURSION 0);	--in case of recursion depth err
	OPEN @curs;
go


--Создать хранимую процедуру, вызывающую процедуру п.1., 
--осуществляющую прокрутку возвращаемого курсора и выводящую сообщения, 
--сформированные из записей при выполнении условия, заданного еще одной пользовательской функцией.

if OBJECT_ID(N'dbo.external_proc', N'P') is not null
	drop procedure dbo.external_proc
go
create procedure dbo.external_proc
AS
	declare @ext_curs cursor;
	declare @t_bkname varchar(254);
	declare @t_auth varchar(35);
	declare @t_antiq int;

	exec dbo.sub_proc @curs = @ext_curs OUTPUT;

	--always use fetch_next before actual fetch according to MSDN
	FETCH NEXT FROM @ext_curs INTO @t_bkname, @t_auth, @t_antiq;
	print 'First Fetch: "' + @t_bkname + '"'

	WHILE (@@FETCH_STATUS = 0)
	BEGIN
		IF (dbo.compare_years(@t_antiq, 40) = 1)
			print '"' + @t_bkname + '" was written by ' + @t_auth + ', ' + CAST(@t_antiq as varchar) + ' years ago'
		FETCH NEXT FROM @ext_curs
		INTO @t_bkname, @t_auth, @t_antiq;
	END

	CLOSE @ext_curs;
	DEALLOCATE @ext_curs;
GO

exec dbo.external_proc
go


--Модифицировать хранимую процедуру п.2. таким образом, 
--чтобы выборка формировалась с помощью табличной функции.

IF EXISTS (SELECT * FROM sysobjects WHERE id = object_id(N'dbo.get_classic') 
    AND xtype IN (N'FN', N'IF', N'TF'))
    DROP FUNCTION dbo.get_classic
go
--table inline function: CREATE FUNCTION ... RETURNS TABLE AS RETURN ( SELECT ... );
create function dbo.get_classic()
	--not inline function:
	returns @tt table
	(
		classic_book_name nvarchar(254),
		classic_book_year int
	)
	as
	begin
		insert @tt
			select books.book_name, books.year
			from dbo.books
			where books.year < 1980
		return
	end
go

alter procedure dbo.sub_proc
	@curs cursor VARYING OUTPUT
as
begin
	set nocount on;
	set @curs = cursor
	scroll static for
		select classic_book_name, classic_book_year
		from dbo.get_classic();		--table function
	open @curs;
end
go

declare @another_curs cursor;

EXEC dbo.sub_proc @curs = @another_curs OUTPUT;
fetch next from @another_curs;
while (@@FETCH_STATUS = 0)
begin
	fetch next from @another_curs;
end
close @another_curs;
deallocate @another_curs;
go

