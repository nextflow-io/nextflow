(mail-page)=

# Notifications

This page documents how to handle workflow events and send notifications.

(workflow-handlers)=

## Workflow handlers

(metadata-completion-handler)=

### Completion handler

Due to the asynchronous nature of Nextflow the termination of a script does not correspond to the termination of the running workflow. Thus some information, only available on execution completion, needs to be accessed by using an asynchronous handler.

The `onComplete` event handler is invoked by the framework when the workflow execution is completed. It allows one to access the workflow termination status and other useful information. For example:

```nextflow
workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
```

(metadata-error-handler)=

### Error handler

The `onError` event handler is invoked by Nextflow when a runtime or process error caused the pipeline execution to stop. For example:

```nextflow
workflow.onError {
    println "Error: Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
```

:::{note}
Both the `onError` and `onComplete` handlers are invoked when an error condition is encountered. The first is called as soon as the error is raised, while the second is called just before the pipeline execution is about to terminate. When using the `finish` {ref}`process-error-strategy`, there may be a significant gap between the two, depending on the time required to complete any pending job.
:::

## Mail

The built-in function `sendMail` allows you to send a mail message from a workflow script.

(mail-basic)=

### Basic mail

The mail attributes are specified as named parameters or an equivalent map. For example:

```nextflow
sendMail(
    to: 'you@gmail.com',
    subject: 'Catch up',
    body: 'Hi, how are you!',
    attach: '/some/path/attachment/file.txt'
)
```

which is equivalent to:

```nextflow
mail = [
    to: 'you@gmail.com',
    subject: 'Catch up',
    body: 'Hi, how are you!',
    attach: '/some/path/attachment/file.txt'
]

sendMail(mail)
```

The following parameters can be specified:

`to`
: *Multiple email addresses can be specified separating them with a comma.*
: The mail target recipients.

`cc`
: *Multiple email addresses can be specified separating them with a comma.*
: The mail CC recipients.

`bcc`
: *Multiple email addresses can be specified separating them with a comma.*
: The mail BCC recipients.

`from`
: *Multiple email addresses can be specified separating them with a comma.*
: The mail sender address.

`subject`
: The mail subject.

`charset`
: The mail content charset (default: `UTF-8`).

`text`
: The mail plain text content.

`body`
: The mail body content. It can be either plain text or HTML content.

`type`
: The mail body mime type. If not specified it's automatically detected.

`attach`
: Single file or a list of files to be included as mail attachments.

(mail-advanced)=

### Advanced mail

Another version of `sendMail` allows a more idiomatic syntax:

```nextflow
sendMail {
    to 'you@gmail.com'
    from 'me@gmail.com'
    attach '/some/path/attachment/file.txt'
    attach '/other/path/image.png'
    subject 'Catch up'

    '''
    Hi there,
    Look! Multi-lines
    mail content!
    '''
}
```

The same attributes listed in the table in the previous section are allowed.

:::{tip}
A string expression at the end is implicitly interpreted as the mail body content, therefore the `body` parameter can be omitted as shown above.
:::

:::{tip}
To send an email that includes text and HTML content, use both the `text` and `body` attributes. The first is used for the plain text content, while the second is used for the rich HTML content.
:::

(mail-attachments)=

### Mail attachments

When using the curly brackets syntax, the `attach` parameter can be repeated two or more times to include multiple attachments in the mail message.

Moreover for each attachment it's possible to specify any of the following options:

`contentId`
: Defines the `Content-ID` header field for the attachment.

`disposition`
: Defines the `Content-Disposition` header field for the attachment.

`fileName`
: Defines the `filename` parameter of the `Content-Disposition` header field.

`description`
: Defines the `Content-Description` header field for the attachment.

For example:

```nextflow
sendMail {
    to 'you@dot.com'
    attach '/some/file.txt', fileName: 'manuscript.txt'
    attach '/other/image.png', disposition: 'inline'
    subject 'Sending documents'
    '''
    the mail body
    '''
}
```

(mail-config)=

### Mail configuration

If no mail server configuration is provided, Nextflow tries to send the email by using the external mail command eventually provided by the underlying system (e.g. `sendmail` or `mail`).

If your system does not provide access to none of the above you can configure a SMTP server in the `nextflow.config` file. For example:

```groovy
mail {
    smtp.host = 'your.smtp-server.com'
    smtp.port = 475
    smtp.user = 'my-user'
}
```

See the {ref}`mail scope <config-mail>` section to learn more the mail server configuration options.

### AWS SES configuration

:::{versionadded} 23.06.0-edge
:::

Nextflow supports [AWS SES](https://aws.amazon.com/ses/) native API as an alternative
provider to send emails in place of SMTP server.

To enable this feature add the following environment variable in the launching environment:

```bash
export NXF_ENABLE_AWS_SES=true
```

Make also sure to add the following AWS IAM permission to the AWS user (or role) used to launch the pipeline execution:

```
ses:SendRawEmail
```

## Mail notification

You can use the `sendMail` function with a {ref}`workflow completion handler <metadata-completion-handler>` to notify the completion of a workflow completion. For example:

```nextflow
workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: 'you@gmail.com', subject: 'My pipeline execution', body: msg)
}
```

This is useful to send a custom notification message. Note however that Nextflow includes a built-in notification mechanism which is the most convenient way to notify the completion of a workflow execution in most cases. Read the following section to learn about it.

## Workflow notification

Nextflow includes a built-in workflow notification features that automatically sends a notification message when a workflow execution terminates.

To enable simply specify the `-N` option when launching the pipeline execution. For example:

```bash
nextflow run <pipeline name> -N <recipient address>
```

It will send a notification mail when the execution completes similar to the one shown below:

```{image} _static/workflow-notification-min.png
```

:::{warning}
By default the notification message is sent with the `sendmail` system tool, which is assumed to be available in the environment where Nextflow is running. Make sure it's properly installed and configured. Alternatively, you can provide the SMTP server configuration settings to use the Nextflow built-in mail support, which doesn't require any external system tool.
:::

See the [Mail configuration](#mail-configuration) section to learn about the available mail delivery options and configuration settings.

Read {ref}`Notification scope <config-notification>` section to learn more about the workflow notification configuration details.
