# SQL DB plugin for Nextflow

This plugin provider an extension to implement built-in support for SQL DB access and manipulation in Nextflow scripts. 

It provides the ability to create Nextflow channel from SQL queries and to populate database tables. The current version 
provides out-of-the-box support for the following databases: 

* [H2](https://www.h2database.com)
* [MySQL](https://www.mysql.com/) 
* [MariaDB](https://mariadb.org/)
* [PostgreSQL](https://www.postgresql.org/)
                    
NOTE: THIS IS A PREVIEW TECHNOLOGY, FEATURES AND CONFIGURATION SETTINGS CAN CHANGE IN FUTURE RELEASES.

This repository only holds plugin artefacts. Source code is available at this [link](https://github.com/nextflow-io/nextflow/tree/master/plugins/nf-sqldb).

## Get started 
  
Make sure to have Nextflow 21.08.0 or later. Add the following snippet to your `nextflow.config` file. 

```
plugins {
  id 'nf-sqldb@0.1.0'
}
```
                                                              
The above declaration allow the use of the SQL plugin functionalities in your Nextflow pipelines. See the section 
below to configure the connection properties with a database instance. 

## Configuration

The target database connection coordinates are specified in the `nextflow.config` file using the
`sql.db` scope. The following are available

| Config option 	                    | Description 	                |
|---	                                |---	                        |
| `sql.db.'<DB-NAME>'.url`      | The database connection URL based on Java [JDBC standard](https://docs.oracle.com/javase/tutorial/jdbc/basics/connecting.html#db_connection_url). 
| `sql.db.'<DB-NAME>'.driver`   | The database driver class name (optional).
| `sql.db.'<DB-NAME>'.user`     | The database connection user name.
| `sql.db.'<DB-NAME>'.password` | The database connection password.

For example:

```
sql {
    db {
        foo {
              url = 'jdbc:mysql://localhost:3306/demo'
              user = 'my-user'
              password = 'my-password'
            }
    }
}

```

The above snippet defines SQL DB named *foo* that connects to a MySQL server running locally at port 3306 and
using `demo` schema, with `my-name` and `my-password` as credentials.

## Available operations

This plugin adds to the Nextflow DSL the following extensions that allows performing query and populate database tables.

### fromQuery

The `fromQuery` factory method allows performing a query against a SQL database and creating a Nextflow channel emitting
a tuple for each row in the corresponding result set. For example:

```
ch = channel.sql.fromQuery('select alpha, delta, omega from SAMPLE', db: 'foo')
```

### insertInto

The `insertInto` operator provided by this plugin allows populating a database table with the data emitted
by a Nextflow channels and therefore produced as result by a pipeline process or an upstream operator. For example:

```
channel
    .of('Hello','world!')
    .map( it -> tuple(it, it.length) )
    .sqlInsert( into: 'SAMPLE', columns: 'NAME, LEN', db: 'foo' )

```

The above example creates and performs the following two SQL statements into the database with name `foo` as defined
in the `nextflow.config` file.

```
INSERT INTO SAMPLE (NAME, LEN) VALUES ('HELLO', 5);
INSERT INTO SAMPLE (NAME, LEN) VALUES ('WORLD!', 6);
```

NOTE: the target table (e.g. `SAMPLE` in the above example) must be created ahead.

The following options are available:

| Operator option 	| Description 	                |
|---	            |---	                        |
| `into`            | The database table name into with the data needs to be stored.
| `columns`         | The database table column names to be filled with the channel data. The column names order and cardinality must match the tuple values emitted by the channel. The columns can be specified as a `List` object or a comma-separated value string.
| `statement`       | The SQL `insert` statement to be performed to insert values in the database using `?` as placeholder for the actual values, for example: `insert into SAMPLE(X,Y) values (?,?)`. When provided the `into` and `columsn` parameters are ignored.
| `db`              | The database handle. It must must a `sql.db` name defined in the `nextflow.config` file.


## Query CSV files

The SQL plugin includes the [H2](https://www.h2database.com/html/main.html) database engine that allows the query of CSV files
as DB tables using SQL statements.

For example, create CSV file using the snippet below:

```
cat <<EOF > test.csv
foo,bar
1,hello
2,ciao
3,hola
4,bonjour
EOF
```

To query this file in a Nextflow script use the following snippet:

```nextflow
    channel
          .sql
          .fromQuery("SELECT * FROM CSVREAD('test.csv') where foo>=2;")
          .view()
```


The `CSVREAD` function provided by the H2 database engine allows the access of a CSV file in your computer file system,
you can replace `test.csv` with a CSV file path of your choice. The `foo>=2` condition shows how to define a filtering
clause using the conventional SQL WHERE constrains. 
