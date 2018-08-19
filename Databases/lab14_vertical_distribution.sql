--Задание 14. Создание вертикально фрагментированных таблиц средствами СУБД SQL Server 2012


--1.Создать в базах данных пункта 1 задания 13 таблицы, содержащие вертикально фрагментированные данные.
--2.Создать необходимые элементы базы данных (представления, триггеры), 
	--обеспечивающие работу с данными вертикально фрагментированных таблиц (выборку, вставку, изменение, удаление).



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
 
 
--1.Создать в базах данных пункта 1 задания 13 таблицы, содержащие вертикально фрагментированные данные.
use lab13_1;
go
if OBJECT_ID(N'dbo.books', N'U') is not null
    drop table dbo.books;
go
create table dbo.books (
    id              int,	--здесь не должно быть identity
							--т.к. в триггере (ins) происходит вставка значения identity в 2 отдельные таблицы.
							--если во время вставки связь с одной таблицей потеряется, то в таблицу не будет занесена запись с верным identity.
							--в следующий раз будет выбрано неверное значение identity для этой таблицы
							--в итоге i-я запись в одной таблице будет соответствовать j-й записи в другой
    title           varchar(254),
 
    PRIMARY KEY (id)
    );
go
 
 
use lab13_2;
go
if OBJECT_ID(N'dbo.books', N'U') is not null
    drop table dbo.books;
go
create table dbo.books (
    id              int,
    writer          varchar(35),
 
    PRIMARY KEY (id),
    );
go
 
 
--2.Создать необходимые элементы базы данных (представления, триггеры), 
	--обеспечивающие работу с данными вертикально фрагментированных таблиц (выборку, вставку, изменение, удаление).
use lab13_1;
go
if OBJECT_ID(N'vertical_dist_v', N'V') is not null
    drop view vertical_dist_v;
go
create view vertical_dist_v as
    select one.id, one.title, two.writer
        from lab13_1.dbo.books as one,
            lab13_2.dbo.books as two
        where one.id = two.id
go
 
 
--==============================================================
-------------------------INSERT-ON-VIEW-------------------------
--==============================================================
if OBJECT_ID(N'dbo.dist_ins', N'TR') is not null
    drop trigger dbo.dist_ins
go
create trigger dbo.dist_ins
    on dbo.vertical_dist_v
    instead of insert
    as
    begin
 
        insert into lab13_1.dbo.books(id, title)
            select id, title
                from inserted
 
        insert into lab13_2.dbo.books(id, writer)
            select id, writer
                from inserted
    end
go
 
insert into dbo.vertical_dist_v values
    (1, 'the lord of the rings', 'john tolkien'),
    (2, 'the hobbit, or there and back again', 'john tolkien'),
    (3, 'a dream of spring', 'martin scorsese'),
    (4, 'a game of thrones', 'george martin'),
    (5, 'random book', 'unknown');
 
select * from dbo.vertical_dist_v
select * from lab13_1.dbo.books
select * from lab13_2.dbo.books
 
 
--==============================================================
-------------------------DELETE-ON-VIEW-------------------------
--==============================================================
if OBJECT_ID(N'dbo.dist_del', N'TR') is not null
    drop trigger dbo.dist_del
go
create trigger dbo.dist_del
    on dbo.vertical_dist_v
    instead of delete
    as
    begin
       
        delete lab13_1.dbo.books
            where id in (select d.id
                from deleted as d)
 
        delete lab13_2.dbo.books
            where id in (select d.id
                from deleted as d)
    end
go
 
delete dbo.vertical_dist_v
    where writer = 'unknown'
 
select * from dbo.vertical_dist_v
select * from lab13_1.dbo.books
select * from lab13_2.dbo.books
 
 
--==============================================================
-------------------------UPDATE-ON-VIEW-------------------------
--==============================================================
if OBJECT_ID(N'dbo.dist_upd', N'TR') is not null
    drop trigger dbo.dist_upd
go
create trigger dbo.dist_upd
    on dbo.vertical_dist_v
    instead of update
    as
    begin
 
        if UPDATE(id)
            raiserror('[UPD TRIGGER]: "id" cant be modidfied', 16, 1)
        else
        begin
 
            update lab13_1.dbo.books
                set books.title = (select title from inserted where inserted.id = books.id)
                    where books.id = (select id from inserted where inserted.id = books.id)
 
            update lab13_2.dbo.books
                set books.writer = (select writer from inserted where inserted.id = books.id)
                    where books.id = (select id from inserted where inserted.id = books.id)
        end
    end
go
 
update dbo.vertical_dist_v
    set writer = 'george martin'
    where writer like '%martin%'
 
select * from dbo.vertical_dist_v
select * from lab13_1.dbo.books
select * from lab13_2.dbo.books

