# Хранилища данных

- присоединение pgadmin'а к postgres:

```bash
docker network create -d bridge pgadmin-net

docker run -p 3001:80 \
    -e "PGADMIN_DEFAULT_EMAIL=ics@a.com" \
    -e "PGADMIN_DEFAULT_PASSWORD=ics" \
    --network homyak_pgadmin-net 
    dpage/pgadmin4

```

- проверка postgres

```postgres-sql
create table Tests
(
    id int primary key,
    foo varchar(20)
)


insert into Tests values
(1, 'hello'),
(2, 'world')

delete from Tests

select * from Tests
```

- проверка SQL server (creds: `sa/Pa55w0rd`)

```postgres-psql
CREATE DATABASE TestDB

SELECT Name from sys.Databases
GO

USE TestDB
CREATE TABLE Inventory (id INT, name NVARCHAR(50), quantity INT)
INSERT INTO Inventory VALUES (1, 'banana', 150);
INSERT INTO Inventory VALUES (2, 'orange', 154);
GO

select * from Inventory
```


- oracle (`jdbc:oracle:thin:@localhost:1521:xe`, `system/oracle`):

```oracle-sql
select banner from v$version where rownum = 1;

ALTER SESSION SET CURRENT_SCHEMA = "SYSTEM";

CREATE TABLE customers
( customer_id number(10) NOT NULL,
  customer_name varchar2(50) NOT NULL,
  city varchar2(50)
);

insert into customers
(customer_id, customer_name, city) VALUES
(1, 'anthony', 'moscow');

select * from customers;
```
