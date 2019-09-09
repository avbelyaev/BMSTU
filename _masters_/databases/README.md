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
