pipeline {

    agent any
    tools {
        maven 'mvn-default'
        jdk   'openjdk-8'
    }

    options {
        timeout(time: 30, unit: 'MINUTES')
    }

    stages {

        stage('Cleaning') {
            steps {
                sh 'git clean -fdx'
                sh 'mvn clean'
            }
        }

        stage('Build') {
            steps {
                sh 'mvn package test install checkstyle:checkstyle'
                checkstyle pattern: 'hipparchus-*/target/checkstyle-result.xml'
                junit 'hipparchus-*/target/surefire-reports/*.xml'
                jacoco execPattern:'**/jacoco.exec', classPattern: '**/classes', sourcePattern: '**/src/main/java', exclusionPattern: 'hipparchus-migration/*'
            }
        }

    }

    post {
        always {
            archiveArtifacts artifacts: 'hipparchus-*/target/*.jar', fingerprint: true
            script {
                if ( env.BRANCH_NAME ==~ /^release-[.0-9]+$/ ) {
                    archiveArtifacts artifacts: 'target/*.zip', fingerprint: true
                }
            }
        }
    }
}
