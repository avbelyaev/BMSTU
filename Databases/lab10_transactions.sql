--Лабораторная работа №10. Режимы выполнения транзакций


--1.Исследовать и проиллюстрировать на примерах различные уровни изоляции транзакций MS SQL Server, 
	--устанавливаемые с использованием инструкции SET TRANSACTION ISOLATION LEVEL.



--сначала запускается одна из транзакций QUERY1, во время ее выполнения параллельно запускается одна из транзакций QUERY2.


QUERY1:
if OBJECT_ID('dbo.cardholders') is not null
    drop table dbo.cardholders
go
create table dbo.cardholders 
(
	id			int,
	name		varchar(35),
	cardtype	varchar(35),
	balance		int,

	PRIMARY KEY (id)
);
go

insert into dbo.cardholders values
	(1, 'George', 'Visa', 1000),
	(2, 'Steve', 'Mastercard, bescenno', 1337),
	(3, 'Samwise the Brave', 'American Express', 1500)
go

--«грязное» чтение(dirty read) – чтение транзакцией записи, измененной другой транзакцией, 
	--при этом эти изменения еще не зафиксированы;

--–невоспроизводимое чтение(non-repeatable read) – при повторном чтении транзакция обнаруживает 
	--измененные или удаленные данные, зафиксированные другой завершенной транзакцией; (update/delete)

--–фантомное чтение(phantom read) – при повторном чтении транзакция обнаруживает новые строки, 
	--вставленные другой завершенной транзакцией; (insert)


--============================================================
--READ UNCOMMITED
--Черновое чтение (Dirty read)
--читатели могт считывать данные незваршенной транзакции
--процесса-писателя
--============================================================

begin transaction
	select * from dbo.cardholders
	update dbo.cardholders 
		set balance = 90 
		where id = 1
	waitfor delay '00:00:05'
	select * from dbo.cardholders
commit




--============================================================
--READ COMMITED
--Подтвержденное чтение
--читатели не могут считывать данные незавершенной транзакции, 
--но писатели могут изменять уже прочитанные данные.
--Если таблица захвачена, то прочитать данные 
--можно только полсле коммита
-->предотвращает грязное чтение
--============================================================

begin transaction
	select * from dbo.cardholders
	update dbo.cardholders 
		set balance = 90 
		where id = 1
	waitfor delay '00:00:05'
		--моментальный результат у читателя, т.к. таблица не захвачена
		--только upd/del захватывают таблицу
	
	select * from dbo.cardholders
commit




--============================================================
--REPEATABLE READ
--Повторяемое чтение
--повторное чтение данных вернет те же значения, 
--что были и в начале транзакции. 
--При этом писатели могут вставлять новые записи, 
--имеющие статус фантома при незавершенной транзакции.
-->предотвращает невоспроизводимое чтение
--============================================================

set transaction isolation level 
	repeatable read
begin transaction
	select * from dbo.cardholders where id in (1, 2)
	waitfor delay '00:00:05'
	select * from dbo.cardholders where id in (1, 2)
		--ins разрешен
		--изменения, затрагивающие данные, выбранные в транзакции, запрещены
		--если данные (строки) не захвачены, их можно менять
	rollback --?
select * from dbo.cardholders




--============================================================
--SERIALIZABLE
--Сериализуемость
--максимальный уровень изоляции, 
--гарантирует неизменяемость данных другими процессами 
--до завершения транзакции.
-->предотвращает фантомное чтение
--============================================================

set transaction isolation level 
	serializable
begin transaction
	select * from dbo.cardholders where id in (1,2)
	waitfor delay '00:00:05'
	select * from dbo.cardholders where id in (1,2)
		--только строки в диапазоне (1 .. 2) PRIMARY KEY будут захвачены
		--остальные строки таблицы можно изменять
	rollback --?
select * from dbo.cardholders





-------------------------------------------------------------------------
QUERY2:


--READ UNCOMMITED

set transaction isolation level 
	read uncommitted
select * from dbo.cardholders



--READ COMMITED

set transaction isolation level 
	read committed
select * from dbo.cardholders



--REPEATABLE READ

update dbo.cardholders set balance = 90 where id = 3
insert into dbo.cardholders values (9, 'John', 'Visa', 0)



--SERIALIZABLE

update dbo.cardholders set balance = 90 where id = 2

