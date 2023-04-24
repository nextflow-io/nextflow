(mail-page)=

# Mail & Notifications

## Mail message

The built-in function `sendMail` allows you to send a mail message from a workflow script.

(mail-basic)=

### Basic mail

The mail attributes are specified as named parameters or providing an equivalent associative array as argument. For example:

```groovy
sendMail(
    to: 'you@gmail.com',
    subject: 'Catch up',
    body: 'Hi, how are you!',
    attach: '/some/path/attachment/file.txt'
)
```

which is equivalent to:

```groovy
mail = [
    to: 'you@gmail.com',
    subject: 'Catch up',
    body: 'Hi, how are you!',
    attach: '/some/path/attachment/file.txt'
]

sendMail(mail)
```

The following parameters can be specified:

| Name          | Description                                                            |
| ------------- | ---------------------------------------------------------------------- |
| to {sup}`*`   | The mail target recipients.                                            |
| cc {sup}`*`   | The mail CC recipients.                                                |
| bcc {sup}`*`  | The mail BCC recipients.                                               |
| from {sup}`*` | The mail sender address.                                               |
| subject       | The mail subject.                                                      |
| charset       | The mail content charset (default: `UTF-8`).                           |
| text          | The mail plain text content.                                           |
| body          | The mail body content. It can be either plain text or HTML content.    |
| type          | The mail body mime type. If not specified it's automatically detected. |
| attach        | Single file or a list of files to be included as mail attachments.     |

`*` Multiple email addresses can be specified separating them with a comma.

(mail-advanced)=

### Advanced mail

Another version of `sendMail` allows a more idiomatic syntax:

```groovy
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

Moreover for each attachment it's possible to specify one or more of the following optional attributes:

| Name        | Description                                                                 |
| ----------- | --------------------------------------------------------------------------- |
| contentId   | Defines the `Content-ID` header field for the attachment.                   |
| disposition | Defines the `Content-Disposition` header field for the attachment.          |
| fileName    | Defines the `filename` parameter of the "Content-Disposition" header field. |
| description | Defines the `Content-Description` header field for the attachment.          |

For example:

```groovy
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

## Mail notification

You can use the `sendMail` function with a {ref}`workflow completion handler <metadata-completion-handler>` to notify the completion of a workflow completion. For example:

```groovy
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

```{image} images/workflow-notification-min.png
```

:::{warning}
By default the notification message is sent with the `sendmail` system tool, which is assumed to be available in the environment where Nextflow is running. Make sure it's properly installed and configured. Alternatively, you can provide the SMTP server configuration settings to use the Nextflow built-in mail support, which doesn't require any external system tool.
:::

See the [Mail configuration](#mail-configuration) section to learn about the available mail delivery options and configuration settings.

Read {ref}`Notification scope <config-notification>` section to learn more about the workflow notification configuration details.
