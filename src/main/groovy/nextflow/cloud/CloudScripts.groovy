package nextflow.cloud

class CloudScripts {
    static String scriptCreateUser(String userName, String key) {
        """\
        useradd -m -s /bin/bash $userName
        mkdir ~$userName/.ssh
        echo "${key.trim()}" > ~$userName/.ssh/authorized_keys
        chmod 700 ~$userName/.ssh
        chmod 600 ~$userName/.ssh/authorized_keys
        chown -R $userName:$userName ~$userName/.ssh
        egrep -i "^wheel:" /etc/group > /dev/null && usermod -aG wheel $userName
        egrep -i "^docker:" /etc/group > /dev/null && usermod -aG docker $userName
        chmod +x /etc/sudoers
        echo '${userName} ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
        chmod -x /etc/sudoers
        """
                .stripIndent()
    }
}
