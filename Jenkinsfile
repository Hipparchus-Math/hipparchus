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
                jacoco execPattern:'**/jacoco.exec', classPattern: 'hipparchus-clustering/target/classes,hipparchus-core/target/classes,hipparchus-fft/target/classes,hipparchus-filtering/target/classes,hipparchus-fitting/target/classes,hipparchus-geometry/target/classes,hipparchus-ode/target/classes,hipparchus-optim/target/classes,hipparchus-stat/target/classes', sourcePattern: 'hipparchus-clustering/src/main/java,hipparchus-core/src/main/java,hipparchus-fft/src/main/java,hipparchus-filtering/src/main/java,hipparchus-fitting/src/main/java,hipparchus-geometry/src/main/java,hipparchus-ode/src/main/java,hipparchus-optim/src/main/java,hipparchus-stat/src/main/java'
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
