
## Lab 1

### 1. Сервер

1. Запустить сервер.
`docker-compose up`

2. Остановить сервер
Ctrl-C

3. Посмотреть журнал.

журнал - в `./postgres-logs`

пример записей:
```
2019-10-08 23:07:19.192 UTC [24] LOG:  database system was not properly shut down; automatic recovery in progress
2019-10-08 23:07:19.289 UTC [24] LOG:  invalid record length at 0/1638080: wanted 24, got 0
2019-10-08 23:07:19.289 UTC [24] LOG:  redo is not required
2019-10-08 23:07:19.508 UTC [1] LOG:  database system is ready to accept connections
2019-10-08 23:16:46.606 UTC [61] ERROR:  database "ics" already exists
2019-10-08 23:16:46.606 UTC [61] STATEMENT:  CREATE DATABASE ics
        WITH 
        OWNER = ics
        ENCODING = 'UTF8'
        LC_COLLATE = 'en_US.utf8'
        LC_CTYPE = 'en_US.utf8'
        TABLESPACE = myts1
        CONNECTION LIMIT = -1;
2019-10-08 23:17:26.265 UTC [61] ERROR:  cannot drop the currently open database
2019-10-08 23:17:26.265 UTC [61] STATEMENT:  DROP DATABASE ics;
```

4. Определить версию сервера

```postgres-sql
SELECT verison();
-- PostgreSQL 10.10 (Debian 10.10-1.pgdg90+1) 
--   on x86_64-pc-linux-gnu, compiled by gcc
```




### 2. Настройки сервера

внутренности прокинуты наружу в `./postgres-data` -> в него можно не заходить! 

1. Посмотрите файл настройки сервера

в контейнере под адресу `/var/lib/postgresql/data/postgresql.conf`. Разделы конфига:
- аутнетификация и авторзация
- файлы
- WAL
- соединение
- оптимизация 
- обработка ошибок
- vacuum

2. Посмотреть содержание pg_log, pg_сlog, pg_хlog

в контейнере под адресу `/var/lib/postgresql/data`

- pg_log - является по умолчанию местом, где хранятся журналы деятельности. 
Они включают в себя сообщения об ошибках, записи о запросах, 
и сообщения во время старта\выключения СУБД.

- pg_xlog (PG10 -> pg_wal) — это место, где PostgreSQL хранит журнал транзакций. 
Этот набор бинарных файлов, с названиями вида '00000001000000000000008E', 
которые содержат образы данных последних транзакций. 

- pg_clog (PG10 -> pg_xact) - содержит журналы метаданных транзакций. 
Этот журнал говорит серверу, какие транзакции завершены, а какие нет. 
Если вы когда-нибудь удалите файлы из pg_clog, вы можете смело удалить и весь каталог 
базы данных. Не существует способа восстановить базу данных без этих журналов.

3. Посмотреть настройки конфигурации 

посмотрели

4. Посмотреть значение параметра shared_buffers

`shared_buffers = 128MB			# min 128kB`

5. Увеличите в три раза 'heap memory per user' процедурой sp_configure

какой нафиг sp_configure. это вообще SQL server

6. Изменить параметр файла конфигурации

Изменили `shared_buffers = 256MB` в `custom-pg.conf`. Рестартанули контейнер



### FAQ

- в PG `use schema1` == `SET search_path TO schema1`