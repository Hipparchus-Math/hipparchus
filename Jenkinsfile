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
            }
        }

        stage('Build') {
            steps {
                sh 'mvn -D maven.test.skip=true package install checkstyle:checkstyle'
                checkstyle pattern: 'hipparchus-*/target/checkstyle-result.xml'
            }
        }

        stage('Test') {
            steps {
                sh 'mvn test'
                junit 'hipparchus-*/target/surefire-reports/*.xml'
                jacoco execPattern:'hipparchus-*/target/**.exec', classPattern: '**/classes', sourcePattern: '**/src/main/java'
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
