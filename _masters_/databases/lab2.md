
### 1. Сервер

#### 1. Запустить сервер.
`docker-compose up`

#### 2. Остановить сервер
Ctrl-C

#### 3. Посмотреть журнал.
журнал - в `./postgres-logs`

#### 4. Определить версию сервера
```postgres-sql
SELECT verison();
-- PostgreSQL 10.10 (Debian 10.10-1.pgdg90+1) on x86_64-pc-linux-gnu, compiled by gcc
```




### 2. Настройки сервера

#### 1. Посмотрите файл настройки сервера
в контейнере под адресу `/var/lib/postgresql/data/postgresql.conf`

#### 2. Посмотреть содержание pg_log, pg_сlog, pg_хlog

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

#### 3. Посмотреть настройки конфигурации 
посмотрели

#### 4. Посмотреть значение параметра shared_buffers

`shared_buffers = 128MB			# min 128kB`

#### 5. Увеличите в три раза 'heap memory per user' процедурой sp_configure

какой нафиг sp_configure. это вообще SQL server

#### 6. Изменить параметр файла конфигурации

Изменили `shared_buffers = 256MB` в `custom-pg.conf`. Рестартанули контейнер




### 3. Создать tablespace

#### 1. Создать tablespace “myts1” (myts2) в директории MyDB1 (MyDB2) с пользователем postgres

- в контейнере:
 
```bash
mkdir /home/space1 && chown postgres /home/space1
mkdir /home/space2 && chown postgres /home/space2
```

- по отдельности запускаем:

```postgres-sql
CREATE TABLESPACE myts1 OWNER ics LOCATION '/home/space1';
CREATE TABLESPACE myts2 OWNER ics LOCATION '/home/space2';
```



CREATE TABLESPACE space1 OWNER ics LOCATION '/home/space1';

-- заходим в контейнер `docker exec -it <container-id> bash
-- проверяем, что в /home/space1 куча мусора


