/*
 * Licensed to the Hipparchus project under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The Hipparchus project licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.hipparchus.complex;

import java.util.Comparator;

/**
 * Comparator for Complex Numbers.
 * 
 */
public class ComplexComparator implements Comparator<Complex> {

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	public int compare(Complex o1, Complex o2) {
		if (o1 == null && o2 != null) {
			return -1;
		} else if (o1 != null && o2 == null) {
			return 1;
		} else if (o1 == null && o2 == null) {
			return 0;
		}

		if (o1.getReal() < o2.getReal()) {
			return -1;
		} else {
			if (o1.getReal() == o2.getReal()) {
				if (o1.getImaginary() == o2.getImaginary()) {
					return 0;
				} else if (o1.getImaginary() < o2.getImaginary()) {
					return -1;
				} else {
					return 1;
				}
			} else {
				return 1;
			}
		}
	}

}
