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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class ComplexComparatorTest {

    private final ComplexComparator comp = new ComplexComparator();
    private Complex o1 = new Complex(1, 1);
    private Complex o11 = new Complex(1, 1);
    private Complex o2 = new Complex(1, 0);
    private Complex o3 = new Complex(2, 0);

    @Test
    public void test() {
        assertEquals("ok", comp.compare(o1, o2), -1 * comp.compare(o2, o1));
    }

    @Test
    public void test2() {
        assertEquals("ok",
                     ((comp.compare(o1, o2) > 0) && (comp.compare(o2, o3) > 0)),
                     comp.compare(o1, o3) > 0);
    }

    @Test
    public void test3() {
        assertEquals("ok", ((comp.compare(o1, o11) == 0)),
                     comp.compare(o1, o2) == comp.compare(o11, o2));
    }

    @Test
    public void test4() {
        assertTrue("ok", (comp.compare(o1, o11) == 0));
    }
}
