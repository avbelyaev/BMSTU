-- 20. Управление журналом транзакций
-- Посмотрите журнал транзакций
-- Создайте процедуру последнего порога для журнала транзакций

SELECT * FROM information_schema.role_table_grants;


-- 21. Авторизация
-- https://www.digitalocean.com/community/tutorials/how-to-use-roles-and-manage-grant-permissions-in-postgresql-on-a-vps--2

-- 1. Задать новую роль (My_role), которая наследует роль sso_role и без пароля.
-- Для этой роли для БД Apt_Ware разрешить
--   - пользоваться таблицей и
--   - разрешить создавать таблицы в этой базе данных
CREATE ROLE sso_role;

CREATE ROLE my_role LOGIN INHERIT;


-- 2. Создать для себя учетную запись Student назначить роль My_role.
-- 3. Создать в BD MyTest группу MyFrands в MyTest
