pipeline {

    agent any
    tools {
        maven 'mvn-default'
        hdk   'openjdk-8'
    }

    options {
        timeout(time: 30, unit: 'MINUTES')
    }

    stages {

        stage('Cleaning') {
            steps {
                sh 'git clean -fdx'
            }
        }

        stage('Build') {
            steps {
                sh 'mvn compile'
                checkstyle pattern: 'target/checkstyle-result.xml'
            }
        }

        stage('Test') {
            steps {
                sh 'mvn test'
                junit 'target/surefire-reports/*.xml'
                jacoco execPattern:'target/**.exec', classPattern: '**/classes', sourcePattern: '**/src/main/java'
            }
        }

        stage('Deploy') {
            steps {
                sh 'mvn install'
            }
        }
    }

    post {
        always {
            archiveArtifacts artifacts: 'target/*.jar', fingerprint: true
            script {
                if ( env.BRANCH_NAME ==~ /^release-[.0-9]+$/ ) {
                    archiveArtifacts artifacts: 'target/*.zip', fingerprint: true
                }
            }
        }
    }
}
