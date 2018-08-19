--Лабораторная работа №12. Программирование клиентского приложения доступа к данным


--1.Создать пользовательское приложение с использованием технологии ADO.NET для выполнения следующих операций:
	--•просмотра содержимого таблиц;
	--•вставки, удаления, изменения данных в этих таблицах.
--2.Сравнить на практике использование для некоторой таблицы:
	--•несвязного уровня ADO.NET;
	--•связного уровня ADO.NET.
--3.Для хранения строки подключения к источнику данных использовать конфигурационный файл приложения (app.config).



use master;
go
if DB_ID (N'lab12ado') is null
    create database lab12ado
        on (
            NAME = lab12adodat,
            FILENAME = 'C:\Users\anthony\Documents\DB_Labs\lab12ado\lab12adodat.mdf',
            SIZE = 5, MAXSIZE = UNLIMITED, FILEGROWTH = 5 )
        log on (
            NAME = lab12adolog,
            FILENAME = 'C:\Users\anthony\Documents\DB_Labs\lab12ado\lab12adolog.ldf',
            SIZE = 5, MAXSIZE = 20, FILEGROWTH = 5 );
go
 

use lab12ado;
go
if OBJECT_ID(N'dbo.test', N'U') is not null
    drop table dbo.test
go
create table dbo.test (
    id			int PRIMARY KEY,
    name		varchar(35),
	email		varchar(254)
    );
go

if OBJECT_ID(N'dbo.books', N'U') is not null
    drop table dbo.books;
go
create table dbo.books (
    title           varchar(254),
    year            int,
    rating          int
        CONSTRAINT DF_books_rating DEFAULT(0),
    );
go

insert into dbo.test values
	(1337, 'anthony', 'anthony@mail.address'),
	(2, 'me', 'me@mail.ru'),
	(3, 'george martin', 'gmartin@7k.com'),
	(1007, 'Q', 'Q@SIS.co.uk'),
	(6000, 'sumail', 'su@mail'),
	(8800, '555', '35@35');

insert dbo.books values
    ('the lord of the rings', 1954, 10),
    ('the hobbit, or there and back again', 1937, 8),
    ('a game of thrones', 1996, 7),
    ('a dance with dragons', 2011, 9),
    ('1984', 1949, 6),
    ('the martian', 2011, 8),
    ('fahrenheit 451', 1953, 6),
    ('the war of the worlds', 1897, 7),
    ('the da vinci code', 2003, 6); 






//app.config
<?xml version="1.0" encoding="utf-8" ?>
<configuration>
    <startup> 
        <supportedRuntime version="v4.0" sku=".NETFramework,Version=v4.5" />
    </startup>

    <connectionStrings>
      <add name="AdoConnString"
           providerName="System.Data.SqlClient"
           connectionString="Server = (localdb)\u5; 
                Initial Catalog = Lab12ado; 
                Integrated Security = true"/>
    </connectionStrings>

</configuration>



//Program.cs
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Data.SqlClient;
using System.Configuration;
using System.Data;

namespace ADOnet
{
    class Program
{

        static int check_tablename(string table, SqlConnection connection) {
			connection.Open();
            DataTable dTable = connection.GetSchema("TABLES", new string[] { null, null, table });

            if (dTable.Rows.Count > 0) {
				connection.Close();
                return 1;
			}
            Console.WriteLine("table does not exist");
			connection.Close();
            return -1;
        }


		static int check_colname(string colname, string table, SqlConnection connection)
		{
			connection.Open();
			DataTable dTable = connection.GetSchema("TABLES", new string[] { null, null, table });

			if (dTable.Rows.Count > 0) {

				string q = string.Format("select * from {0}", table);
                SqlCommand command = new SqlCommand(q, connection);
                SqlDataReader reader = command.ExecuteReader();
				if (reader.Read()) {

					int fieldCount = reader.FieldCount;
					string[] fields = new string[fieldCount];

					for (int j = 0; j < fieldCount; j++) {
						fields[j] = reader.GetName(j);
					}

					if (fields.Contains<String>(colname)) {
						reader.Close();
						connection.Close();
						return 1;
					}
				}
				reader.Close();
			}
			
			Console.WriteLine("column does not exist");
			connection.Close();
			return -1;
		}
		

        static int run(SqlCommand command, SqlConnection connection) {
            try {
                command.ExecuteNonQuery();

            } catch (SqlException e) { 
                Console.WriteLine(e.Message);
                return -1;
            }
			return 0;
        }


        static void exec(SqlCommand command) {
            try {

				SqlDataReader reader = command.ExecuteReader();
				do {
					while (reader.Read()) {
						Object[] values = new Object[reader.FieldCount];
						int fieldCount = reader.GetValues(values);

						for (int i = 0; i < fieldCount; i++) {
							Console.Write("\t\t{0}", values[i]);
						}
						Console.WriteLine();
					}

				} while (reader.NextResult());

				reader.Close();


            } catch (Exception e) { Console.WriteLine(e.Message); }

            Console.WriteLine();
        }


        static void PrintTable(DataTable dt) {
            DataTableReader dtReader = dt.CreateDataReader();

            while (dtReader.Read()) {
                for (int i = 0; i < dtReader.FieldCount; i++) {
                    Console.Write("\t{0}", dtReader.GetValue(i).ToString());
                }
                Console.WriteLine();
            }
            dtReader.Close();
        }

//---------------------------------------------------------
//--------------------------MAIN---------------------------
//---------------------------------------------------------


		static void Main(string[] args)
		{
			try 
			{



				SqlConnection connection;

				try {
					string connectionString = ConfigurationManager.ConnectionStrings["AdoConnString"].ConnectionString;
					connection = new SqlConnection(connectionString);
				}
				catch (Exception e) { Console.WriteLine(e.Message); Console.ReadLine(); return; }

				Console.WriteLine("Connection has been opened.");
				Console.WriteLine("Cвязный уровень");
				Console.WriteLine("type 'help' for help");

				string com;

				while (true)
				{
					string table, column, q;
					Console.Write(">");
					com = Console.ReadLine();

					if (com.Equals("")) 
						continue;
//=========================================================
//------------------------SELECT---------------------------
//=========================================================
					if ('s' == com[0])
					{
						Console.Write("->tables | * | column:");
						string subcom = Console.ReadLine();

						if (subcom.Equals("")) 
							continue;

						if ('t' == subcom[0])
						{
							q = "select * from INFORMATION_SCHEMA.TABLES";
							connection.Open();
							SqlCommand cmd = new SqlCommand(q, connection);

							exec(cmd);
							connection.Close();
						}

						if ('*' == subcom[0])
						{
							Console.Write("->tablename:");
							table = Console.ReadLine();
							if (-1 == check_tablename(table, connection))
								continue;

							connection.Open();
							q = string.Format("select * from {0}", table);
							SqlCommand cmd = new SqlCommand(q, connection);

							exec(cmd);
							connection.Close();
						}

						if ('c' == subcom[0])
						{
							Console.Write("->tablename:");
							table = Console.ReadLine();
							if (-1 == check_tablename(table, connection))
								continue;

							Console.Write("-->column name:");
							column = Console.ReadLine();
							if (-1 == check_colname(column, table, connection))
								continue;

							connection.Open();
							q = string.Format("select {0} from {1}", column, table);
							SqlCommand cmd = new SqlCommand(q, connection);

							exec(cmd);
							connection.Close();
						}

						if (subcom.Equals("q"))
							continue;
					}


//=========================================================
//------------------------INSERT---------------------------
//=========================================================
					if ('i' == com[0])
					{
						Console.Write("->tablename:");
						table = Console.ReadLine();
						if (-1 == check_tablename(table, connection))
							continue;

						connection.Open();
						q = string.Format("select * from {0}", table);
						SqlCommand command = new SqlCommand(q, connection);
						SqlDataReader reader = command.ExecuteReader();

						//get column names as result of simple <select*> query
						if (reader.Read())
						{

							int fieldCount = reader.FieldCount;

							//form query string
							q = string.Format("insert into {0} (", table);
							for (int j = 0; j < fieldCount; j++)
							{
								q += reader.GetName(j) + ", ";
							}
							q = q.Remove(q.Length - 2);
							q += ") values (";
							for (int j = 0; j < fieldCount; j++)
							{
								q += "@" + reader.GetName(j) + ", ";
							}
							q = q.Remove(q.Length - 2);
							q += ")";
							Console.WriteLine("Q:" + q);


							SqlCommand cmd = new SqlCommand(q, connection);

							//get parameters and add to command
							for (int j = 0; j < fieldCount; j++) {
								Console.Write("-->{0}[{1}]:", reader.GetName(j), reader.GetDataTypeName(j));
								Object value = Console.ReadLine();

								cmd.Parameters.AddWithValue("@" + reader.GetName(j), value);
							}

							reader.Close();

							if (0 == run(cmd, connection))
								Console.WriteLine("entry has been inserted");

						} else {
							reader.Close();
						}
						connection.Close();
					}


//=========================================================
//------------------------DELETE---------------------------
//=========================================================
					if ('d' == com[0])
					{
						Console.Write("->tablename:");
						table = Console.ReadLine();
						if (-1 == check_tablename(table, connection))
							continue;

						connection.Open();
						q = string.Format("select * from {0}", table);
						SqlCommand command = new SqlCommand(q, connection);
						SqlDataReader reader = command.ExecuteReader();

						if (reader.Read())
						{

							Object[] del_values = new Object[reader.FieldCount];//array of values to delete
							int[] field_indices = new int[reader.FieldCount];//array of their indices in table
							int delete_col_num = 0, fieldCount = reader.FieldCount;

							q = string.Format("delete {0} where ", table);

							for (int j = 0; j < fieldCount; j++)
							{
								Console.Write("-->del by {0}? Y/N:", reader.GetName(j));
								string t_com = Console.ReadLine();

								if ('n' == t_com[0] || 'N' == t_com[0]) {
									continue;
								}
								else {
									Console.Write("-->{0}[{1}] value:", reader.GetName(j), reader.GetDataTypeName(j));
									del_values[delete_col_num] = Console.ReadLine();

									field_indices[delete_col_num] = j;
									q += reader.GetName(j) + " = @" + reader.GetName(j) + " and ";
									delete_col_num++;
								}
							}

							if (0 == delete_col_num)
							{
								Console.WriteLine("delete parameters are missing");
								reader.Close();
								connection.Close();
								continue;
							}
							q = q.Remove(q.Length - 5);
							Console.WriteLine("Q:" + q);



							SqlCommand cmd = new SqlCommand(q, connection);

							for (int k = 0; k < delete_col_num; k++) {
								cmd.Parameters.AddWithValue("@" + reader.GetName(field_indices[k]), del_values[k]);
							}
							reader.Close();

							if (0 == run(cmd, connection))
								Console.WriteLine("entry(-ies) has been deleted");

						} else {
							reader.Close();
						}
						connection.Close();
					}


//=========================================================
//------------------------UPDATE---------------------------
//=========================================================
					if ('u' == com[0])
					{
						Console.Write("->tablename:");
						table = Console.ReadLine();
						if (-1 == check_tablename(table, connection))
							continue;

						connection.Open();
						q = string.Format("select * from {0}", table);
						SqlCommand command = new SqlCommand(q, connection);
						SqlDataReader reader = command.ExecuteReader();

						if (reader.Read())
						{
							Object[] old_values = new Object[reader.FieldCount];
							int[] old_field_indices = new int[reader.FieldCount];//array of their indices in table
							Object[] new_values = new Object[reader.FieldCount];
							int[] new_field_indices = new int[reader.FieldCount];//array of their indices in table
							int new_col_num = 0, old_col_num = 0, fieldCount = reader.FieldCount;

							q = string.Format("update {0} set ", table);

							for (int j = 0; j < fieldCount; j++) {
								Console.Write("-->modify {0}? Y/N:", reader.GetName(j));
								string t_com = Console.ReadLine();

								if ('n' == t_com[0] || 'N' == t_com[0]) {
									continue;
								}
								else {
									Console.Write("-->{0}[{1}] value:", reader.GetName(j), reader.GetDataTypeName(j));
									new_values[new_col_num] = Console.ReadLine();
									new_field_indices[new_col_num] = j;
									q += reader.GetName(j) + " = @new" + reader.GetName(j) + ", ";
									new_col_num++;
								}
							}
							q = q.Remove(q.Length - 2);

							q += " where ";
							for (int j = 0; j < fieldCount; j++) {
								Console.Write("-->by {0}? Y/N:", reader.GetName(j));
								string t_com = Console.ReadLine();

								if ('n' == t_com[0] || 'N' == t_com[0]) {
									continue;
								} 
								else {
									Console.Write("-->{0}[{1}] value:", reader.GetName(j), reader.GetDataTypeName(j));
									old_values[old_col_num] = Console.ReadLine();
									old_field_indices[old_col_num] = j;
									q += reader.GetName(j) + " = @old" + reader.GetName(j) + " and ";
									old_col_num++;
								}
							}
							q = q.Remove(q.Length - 5);
							Console.WriteLine("Q:" + q);

							if (0 == new_col_num || 0 == old_col_num)
							{
								Console.WriteLine("delete or update parameters are missing");
								reader.Close();
								connection.Close();
								continue;
							}



							SqlCommand cmd = new SqlCommand(q, connection);

							for (int k = 0; k < new_col_num; k++) {
								cmd.Parameters.AddWithValue("@new" + reader.GetName(new_field_indices[k]), new_values[k]);
							}
							for (int k = 0; k < old_col_num; k++) {
								cmd.Parameters.AddWithValue("@old" + reader.GetName(old_field_indices[k]), old_values[k]);
							}

							reader.Close();

							if (0 == run(cmd, connection))
								Console.WriteLine("entry(-ies) has been updated");

						} else {
							reader.Close();
						}
						connection.Close();
					}


//=========================================================
//------------------------EXPLICIT-------------------------
//=========================================================
					if ('e' == com[0])
					{
						connection.Open();
						
						Console.Write("->");
						q = Console.ReadLine();
						SqlCommand cmd = new SqlCommand(q, connection);

						exec(cmd);

						connection.Close();
					}

					if ('h' == com[0])
					{
						Console.WriteLine("Available commands:");
						Console.Write("\tselect\n\t\ttables\n\t\t*\n\t\t\t<table_name>\n\t\tcolumn\n\t\t\t<table_name>\n\t\t\t\t<column_name>\n");
						Console.Write("\tinsert\n\t\t<table_name>\n\t\t\t<value_string>\n");
						Console.Write("\tdelete\n\t\t<table_name>\n\t\t\t<by-fields>\n");
						Console.Write("\tupdate\n\t\t<table_name>\n\t\t\t<fields to modify>\n\t\t\t\t<by-fields>\n");
						Console.WriteLine("\texplicit\n\t\t<explicit_query>\n");
						Console.WriteLine("\tПерейти на Несвязный уроень / Выйти = [q]");
					}

					if (com.Equals("q")) 
						break;

				}

				Console.Write("Выйти = [q]\nПерейти на Несвзный уровень = [Any_key]\nПерейти?:");
				com = Console.ReadLine();

				if (com.Equals("q") && !com.Equals(""))
				{
					Console.WriteLine("Closing connection. App will shut down soon.");
					if (null != connection && connection.State == ConnectionState.Open)
					{
						connection.Close();
					}
					
					return;
				}
				Console.Write("\n\nНесвязный уровень\n");
				Console.WriteLine("type 'help' for help");





				//несвязный уровень

				DataSet cds = new DataSet("cds");

				SqlDataAdapter adapter = new SqlDataAdapter("select * from test", connection);
				SqlCommandBuilder builder = new SqlCommandBuilder(adapter);

				adapter.Fill(cds, "test");

				DataTable dt = cds.Tables["test"];
				
				//table must have primary key for delete/update queries
				try {
					DataColumn[] primary_key = new DataColumn[1];
					primary_key[0] = cds.Tables["test"].Columns["id"];
					dt.PrimaryKey = primary_key;
				} 
				catch (Exception e) { Console.WriteLine(e.Message); }
				



				while (true)
				{

					Console.Write(">");
					com = Console.ReadLine();

					if (com.Equals("")) 
						continue;
//=========================================================
//------------------------SELECT---------------------------
//=========================================================
					if ('s' == com[0])
					{
						Console.WriteLine("select * from test");
						PrintTable(cds.Tables["test"]);
					}

//=========================================================
//------------------------INSERT---------------------------
//=========================================================
					if ('i' == com[0])
					{
						Console.Write("->id:");
						int tid;
						try
						{
							tid = Convert.ToInt32(Console.ReadLine());
						}
						catch (Exception e) { Console.WriteLine(e.Message); continue; }

						Console.Write("->name:");
						string tname = Console.ReadLine();
						Console.Write("->email:");
						string tmail = Console.ReadLine();

						DataRow row = dt.NewRow();
						row[0] = tid;
						row[1] = tname;
						row[2] = tmail;

						try { dt.Rows.Add(row); }
						catch (Exception e) { Console.WriteLine(e.Message); continue; }

						Console.WriteLine("row has been added");
					}

//=========================================================
//------------------------DELETE---------------------------
//=========================================================
					if ('d' == com[0])
					{
						Console.Write("->by PrimaryKey id:");
						DataRow row = null;

						try
						{
							int tid = Convert.ToInt32(Console.ReadLine());
							row = dt.Rows.Find(tid);
						} 
						catch (Exception e) { Console.WriteLine(e.Message); continue; }

						try { row.Delete(); }
						catch (Exception e) { Console.WriteLine(e.Message); continue; }

						Console.WriteLine("row has been deleted");
					}

//=========================================================
//------------------------UPDATE---------------------------
//=========================================================
					if ('u' == com[0])
					{

						Console.Write("->new email:");
						string newmail = Console.ReadLine();

						int tid, indx;
						DataRow row = null;

						Console.Write("->by id:");
						try
						{
							tid = Convert.ToInt32(Console.ReadLine());
							row = dt.Rows.Find(tid);
							indx = dt.Rows.IndexOf(row);
						}
						catch (Exception e) { Console.WriteLine(e.Message); continue; }

						try { dt.Rows[indx]["email"] = newmail; }
						catch (Exception e) { Console.WriteLine(e.Message); continue; }

						Console.WriteLine("row has been updated");
						continue;
					}

//=========================================================
//---------------------ACCEPT-CHANGES----------------------
//=========================================================
					if ('a' == com[0])
					{
						adapter.Update(cds, dt.TableName);
						Console.WriteLine("changes have been accepted");
					}


					if ('h' == com[0])
					{
						Console.WriteLine("Available commands:");
						Console.Write("\tselect\n\t\t* from 'test'\n");
						Console.Write("\tinsert\n\t\tinto 'test'\n\t\t\t<value_string>\n");
						Console.Write("\tdelete\n\t\t'test' by\n\t\t\t<primary_key>\n");
						Console.Write("\tupdate\n\t\t'test' set 'emai'\n\t\t\t<new_email>\n\t\t\t\tby\n\t\t\t\t\t<primary_key>\n");
						Console.WriteLine("\tAccept changes to table = [a]\n");
						Console.Write("\tExit = [q]\n");
					}

					if (com.Equals("q")) 
						break;
				}

				Console.WriteLine("Closing connection. App will shut down soon.");
				connection.Close();

		} catch(Exception e) {
			Console.WriteLine(e.Message);
			Console.ReadLine();
			return;
		}
	
		}

	}

}

