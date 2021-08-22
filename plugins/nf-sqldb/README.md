# SQL DB plugin for Nextflow

This plugin provider an experimental extension to implement 
built-in support for SQL DB in Nextflow. 

It provides the ability to create Nextflow channel from SQL query. 
The query can be made against any SQL DB supporting the JDBC interface
or even CSV/TSV files.

### Get started 

1. Clone and compile the branch `nf-sqldb`.
2. Compile the project with the command

   ```
   make compile
   ```         
3. Create the following CSV snippet 

    ```
    cat <<EOF > test.csv
    foo,bar
    1,hello
    2,ciao
    3,hola
    4,bonjour
    EOF
    ```

4. Save the following Nextflow snippet into a file:    

    ```
    cat <<EOF > sql.nf
    Channel
      .sql
      .fromQuery("SELECT * FROM CSVREAD('test.csv') where foo>=2;")
      .view ()
    EOF
    ```

5. Run it 

    ```
    ./launch.sh run sql.nf -plugins nf-sqldb
   ```


### MySQL examples 

1. Create `nextflow.config` file with this snippet 

   ```
   plugins {
     id 'nf-sqldb'
   }
   
   dataSources {
       mysql {
         url = 'jdbc:mysql://localhost:3306/demo'
         user = 'demo'
         password = 'demo'
       }
   }
   
   ```
  
2. Insert data into a MySQL table

   ```
   channel
     .of('one','two','three')
     .map { it -> tuple(it, it.size()) }
     .sqlInsert( 'insert into sample (name, count) values (?, ?)', dataSource: 'mysql')
   ```
