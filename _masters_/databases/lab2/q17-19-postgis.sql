-- переключаемся на соседний xpggis контейнер

-- 17. PGGis
-- http://strimas.com/spatial/r-postgis-2/

-- 18. Создать БД

-- 1. Создать БД myGis
DROP DATABASE IF EXISTS mygis;
CREATE DATABASE mygis
    WITH
    OWNER = ics
    ENCODING = 'UTF8'
    TABLESPACE = myts1;

/*
переключаемся на БД mygis
psql --username ics --password ics --dbname mygis
 */

-- make sure extension already exists
CREATE EXTENSION postgis;

-- Таблица Rectangle
--     1. Первичный ключ автоинкрементальный (Serial).
--     2. Поле для геометрии MyRec
drop table if exists rectangles;
create table rectangles
(
    id serial primary key,
    name varchar(255),
    geom geometry
);

-- 19. SQL-запросы
-- Insert
-- 1. Добавьте следующие прямоугольники с описнием
insert into rectangles(name, geom) values
('A', 'POLYGON((0 0, 30 0, 30 20, 0 20, 0 0))'),
('B', 'POLYGON((5 5, 10 5, 10 10, 5 10, 5 5))'),
('C', 'POLYGON((20 -5, 35 -5, 35 5, 20 5, 20 -5))'),
('D', 'POLYGON((0 -15, 5 -15, 5 -10, 0 -10, 0 -15))'),
('E', 'POLYGON((20 -20, 30 -20, 30 -15, 20 -15, 20 -20))');

-- Select
-- Сформировать список всех прямоугольников
select name, ST_AsText(geom) from rectangles;
/*
схематично выше были заданы прямоугольники, располагающиеся таким обраом:

---------------------(30, 20)
| A _______(10, 10)  |
|   |  B  |   _____________35, 5)
|   -------   |      |    |
--------------|-------  C |
              -------------
_____(5,-10)
| D |         _____(30, -15)
-----         | E |
              -----
*/

-- Сформировать список всех пересекающихся
-- Make sure to avoid comparing the a polygon to itself.
-- Also avoid comparing a pair of polygon twice
create view intersections as
SELECT r1.name as geom1, r2.name as geom2, st_intersects(r1.geom, r2.geom)
from rectangles r1, rectangles r2
where r1.name < r2.name
    and st_intersects(r1.geom, r2.geom)
    and not st_contains(r1.geom, r2.geom);

-- убеждаемся, что A и C пересекаются. A и B мы исключили,
-- т.к. пересечением считаем пересечение контуров, а не площадей
select * from intersections;
--  geom1 | geom2 | st_intersects
-- -------+-------+---------------
--  A     | C     | t
-- (1 row)

-- Сформировать список всех прямоугольнико лежащих внутри первого
select r.name
from rectangles r
where r.name != 'A'
    and st_contains((select geom
                    from rectangles r1
                    where r1.name = 'A'
        ), r.geom);

-- внутри A лежит только B
select * from intersections;
--  name
-- ------
--  B
-- (1 row)

-- Delete
-- 1. Удалить все строки с непресекающимися прямоугольниками
delete from rectangles
where name not in (
    select geom1 from intersections
    union
    select geom2 from intersections
);

-- остаются только A и C, т.к. на рисунке выше только они пересекались
select name, ST_AsText(geom) from rectangles;
--  name |               st_astext
-- ------+----------------------------------------
--  A    | POLYGON((0 0,30 0,30 20,0 20,0 0))
--  C    | POLYGON((20 -5,35 -5,35 5,20 5,20 -5))
-- (2 rows)
