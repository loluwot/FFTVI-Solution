import java.io.*;
import java.math.BigInteger;
import java.util.*;
public class Main {
    static BufferedReader br;
    static PrintWriter pw;
    static StringTokenizer st;

    static void exit() {
        pw.close();
        System.exit(0);
    }
    private static final BigInteger TWO_64 = BigInteger.ONE.shiftLeft(64);

    public String asUnsignedDecimalString(long l) {
        BigInteger b = BigInteger.valueOf(l);
        if(b.signum() < 0) {
            b = b.add(TWO_64);
        }
        return b.toString();
    }
    static long readLong() throws IOException {
        return Long.parseLong(next());
    }

    static double readDouble() throws IOException {
        return Double.parseDouble(next());
    }

    static int readInt() throws IOException {
        return Integer.parseInt(next());
    }

    static char readChar() throws IOException {
        return next().charAt(0);
    }

    static String nextLine() throws IOException {
        String s = br.readLine();
        if (s == null) {
            exit();
        }
        st = null;
        return s;
    }

    static String next() throws IOException {
        while (st == null || !st.hasMoreTokens()) {
            st = new StringTokenizer(nextLine().trim());
        }
        return st.nextToken();
    }

//    public static class Complex {
//        double a;
//        double b;
//
//        public Complex() {
//            this.a = 0;
//            this.b = 0;
//        }
//
//        public Complex(double a, double b) {
//            this.a = a;
//            this.b = b;
//        }
//
//        public Complex(double a) {
//            this.a = a;
//            this.b = 0;
//        }
//
//        public Complex multiply(Complex other) {
//            double a = this.a * other.a - this.b * other.b;
//            double b = this.a * other.b + this.b * other.a;
//            if (Math.abs(a) < Math.pow(10, -15)) {
//                a = 0;
//            }
//            if (Math.abs(b) < Math.pow(10, -15)) {
//                b = 0;
//            }
//            return new Complex(a, b);
//        }
//
//        public Complex add(Complex other) {
//            double a = this.a + other.a;
//            double b = this.b + other.b;
//            return new Complex(a, b);
//        }
//
//        public Complex negative() {
//            return new Complex(-1 * this.a, -1 * this.b);
//        }
//
//        public Complex conjugate() {
//            return new Complex(this.a, -1 * this.b);
//        }
//
//        public String toString() {
//            return a + " + " + b + "i";
//        }
//    }

//    public static Complex primitiveRoot(int m) {
//        double a = Math.cos(2 * Math.PI / m);
//        double b = Math.sin(2 * Math.PI / m);
//        return new Complex(a, b);
//    }
//
//    public static Complex iPrimitiveRoot(int m) {
//        double a = Math.cos(-2 * Math.PI / m);
//        double b = Math.sin(-2 * Math.PI / m);
//        return new Complex(a, b);
//    }
//
//    public static Complex[] generateRoots(int m) {
//        Complex temp = new Complex(1, 0);
//        Complex primitive = primitiveRoot(m);
//        Complex[] arr = new Complex[m];
//        for (int i = 0; i < m / 2; i++) {
//            arr[i] = temp;
//            arr[i + m / 2] = temp.negative();
//            temp = temp.multiply(primitive);
//        }
//        return arr;
//    }
    public static long [] pad(long[] arr) {
        int low = (int) (Math.ceil(Math.log(arr.length * 2) / Math.log(2)));
        int high = (int) (Math.pow(2, low));
        long [] n = new long [high];
        for (int i = 0; i < arr.length; i++) {
            n[i] = arr[i];
        }
        return n;
    }

//    public static int pad(long arr) {
//        int low = (int) (Math.log(arr)*2/ Math.log(2));
//        int high = (int) (Math.pow(2, low));
//
//        return high;
//    }

    //    public static Complex [] eval(Complex [] arr, int m, Complex primitive){
//        if (m == 1){
//            Complex [] temp = {new Complex(arr[0].a, arr[0].b)};
//            return temp;
//        }
//        else{
//            Complex [] temp = new Complex[m];
//            Complex aodd [] = new Complex [arr.length/2];
//            Complex aeven [] = new Complex [arr.length/2];
//            for (int i = 0; i < arr.length; i++){
//                if (i % 2 ==0){
//                    aeven[i/2] = arr[i];
//                }
//                else{
//                    aodd[i/2] = arr[i];
//                }
//            }
//            Complex [] ffteven = eval(aeven, m/2, primitive.multiply(primitive));
//            Complex [] fftodd = eval(aodd, m/2, primitive.multiply(primitive));
//            Complex x = new Complex(1, 0);
//            for (int i = 0; i < m/2; i++){
//                temp[i] = ffteven[i].add(fftodd[i].multiply(x));
//                temp[i+m/2] = ffteven[i].add((fftodd[i].multiply(x)).negative());
//                x = x.multiply(primitive);
//            }
//            return temp;
//        }
//    }
//    public static void inPlace(Complex[] arr, int m, boolean inverse) {
//        for (int i = 1, j = 0; i < m; i++) {
//            int bit = m >> 1;
//            for (; (j & bit) > 0; bit >>= 1) {
//                j ^= bit;
//            }
//            j ^= bit;
//
//            if (i < j) {
//                Complex temp = arr[i];
//                arr[i] = arr[j];
//                arr[j] = temp;
//            }
//        }
//        Complex primitive;
//        if (inverse) {
//            primitive = iPrimitiveRoot(2);
//        } else {
//            primitive = primitiveRoot(2);
//        }
//        for (int i = 2; i <= m; i <<= 1) {
//            Complex primitive;
//            if (inverse) {
//                primitive = iPrimitiveRoot(i);
//            } else {
//                primitive = primitiveRoot(i);
//            }
//            for (int j = 0; j < m; j += i) {
//                Complex temp = new Complex(1);
//                for (int k = j; k < j + i / 2; k++) {
//                    Complex aeven = arr[k];
//                    Complex aodd = arr[k + i / 2];
//                    arr[k] = aeven.add(aodd.multiply(temp));
//                    arr[k + i / 2] = aeven.add(aodd.multiply(temp).negative());
//                    temp = temp.multiply(primitive);
//                }
//            }
//        }
//    }

    public static List<Long> uniquePrimeFactors(long n) {
        if (n < 1)
            throw new IllegalArgumentException();
        List<Long> result = new ArrayList<>();
        for (long i = 2, end = sqrt(n); i <= end; i++) {
            if (n % i == 0) {
                result.add(i);
                do n /= i;
                while (n % i == 0);
                end = sqrt(n);
            }
        }
        if (n > 1)
            result.add(n);
        return result;
    }

    public static int sqrt(long x) {
        if (x < 0)
            throw new IllegalArgumentException();
        int y = 0;
        for (int i = 1 << 15; i != 0; i >>>= 1) {
            y |= i;
            if (y > 46340 || y * y > x)
                y ^= i;
        }
        return y;
    }

    public static long[] findN(int len) {
        long max = 3037000500l;
        long arr[] = new long[2];
        BigInteger len1 = BigInteger.valueOf(len);

        long start = (max / len + 1) * len + 1;
        int counter = 0;
        for (BigInteger i = BigInteger.valueOf(start); ; i = i.add(len1)) {
            if (i.isProbablePrime(100)) {
                arr[0] = counter + (max / len) + 1;
                arr[1] = i.longValue();
                return arr;
            }
            counter++;
        }
    }

    public static long[] findN2(int len, long firstN, long lastCount) {
        long arr[] = new long[2];
        BigInteger len1 = BigInteger.valueOf(len);
        int counter = 0;
        firstN += len;
        for (BigInteger i = BigInteger.valueOf(firstN); ; i = i.add(len1)) {
            if (i.isProbablePrime(100)) {
                arr[0] = counter + lastCount + 1;
                arr[1] = i.longValue();
                return arr;
            }
            counter++;
        }
    }

//    public static BigInteger findGenerator(List<Long> primeFactors, long modulus) {
//        BigInteger mod = BigInteger.valueOf(modulus);
////        BigInteger totient = BigInteger.valueOf(modulus-1);
//        while (true) {
//            BigInteger possibleGenerator = BigInteger.valueOf((long) (Math.random() * modulus));
//            boolean valid = true;
//            for (long p : primeFactors) {
//                if (possibleGenerator.modPow(BigInteger.valueOf((modulus - 1) / p), mod).longValue() == 1) {
//                    valid = false;
//                    break;
//                }
//            }
//            if (valid) {
//                return possibleGenerator;
//            }
//        }
//    }

//    public static BigInteger findModPrimitive(BigInteger generator, long mod, long k, boolean inverse) {
////        BigInteger primitive = generator.modPow(BigInteger.valueOf(k), BigInteger.valueOf(mod));
////        if (inverse){
////            return primitive.modInverse(BigInteger.valueOf(mod));
////        }
////        return primitive;
////    }
////    public static BigInteger findPrimitive(BigInteger generator, long mod, BigInteger kn, long len, boolean inverse){
////        return findModPrimitive(generator, mod, kn.divide(BigInteger.valueOf(len)).longValue(), inverse);
////    }

    public static void NTT(long [] arr, int m, long ar [], long [] precompute, boolean inverse){
        for (int i = 1, j = 0; i < m; i++) {
            int bit = m >> 1;
            for (; (j & bit) > 0; bit >>= 1) {
                j ^= bit;
            }
            j ^= bit;

            if (i < j) {
                long temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            }
        }
//        System.out.println(Arrays.toString(arr));
        int counter = 0;
        for (int i = 2; i <= m; i <<= 1) {
            long primitive = precompute[counter];

//            System.out.println(primitive + " major primitive");
            for (int j = 0; j < m; j += i) {
                long temp = 1;
                for (int k = j; k < j + (i >> 1); k++) {
//                    System.out.println(temp + " current primitive for length " + i);
                    long aeven = arr[k];
                    long aodd = arr[k + (i >> 1)];
                    arr[k] = (aeven+(aodd*temp)%ar[1])%ar[1];
                    arr[k + (i >> 1)] = (aeven - (aodd*temp)%ar[1] + ar[1])%ar[1];
                    temp = (temp*primitive)%ar[1];
                }

            }
            counter++;
//            System.out.println("Array after length " + i + " :" + Arrays.toString(arr));
        }
        if (inverse){
            for (int i = 0; i < arr.length; i++){
//                arr[i] = (arr[i]*(BigInteger.valueOf(arr.length).modInverse(BigInteger.valueOf(ar[1]))).longValue())%ar[1];
                arr[i] = (arr[i]*euclid(arr.length, ar[1])[1])%ar[1];

            }
        }

    }
//    public static BigInteger [] euclid(BigInteger a, BigInteger b){
//        if (a.equals(BigInteger.ZERO)){
//            BigInteger ret []  = {b, BigInteger.ZERO, BigInteger.ONE};
//            return ret;
//        }
//
//        BigInteger d [] = euclid(b.add(a).mod(a), a);
////        System.out.println("h");
////        System.out.println(Arrays.toString(d));
//        BigInteger m1 = d[2].subtract(b.divide(a).multiply(d[1]));
//        BigInteger m2 = d[1];
//        BigInteger test [] = {d[0], m1, m2};
//        return test;
//    }
        public static long [] euclid(long a, long b){
            if (a == 0){
                long ret []  = {b, 0, 1};
                return ret;
            }

            long d [] = euclid((b+a)%a, a);
            long m1 = d[2] - (b/a) * d[1];
            long m2 = d[1];
            long test [] = {d[0], m1, m2};
            return test;
        }

    public static BigInteger crt (long a1, long a2, long mod1, long mod2, BigInteger bigMod){
        BigInteger mod1b = BigInteger.valueOf(mod1);
        BigInteger mod2b = BigInteger.valueOf(mod2);

        long d[] = euclid (mod1, mod2);
        BigInteger m1 = BigInteger.valueOf(d[1]);
        BigInteger m2 = BigInteger.valueOf(d[2]);
        BigInteger a1b = BigInteger.valueOf(a1);
        BigInteger a2b = BigInteger.valueOf(a2);
        BigInteger ans = a1b.multiply(m2).multiply(mod2b).add(a2b.multiply(m1).multiply(mod1b));
        return ans.mod(bigMod);

    }
    public static void main(String[] args) throws IOException {
        br = new BufferedReader(new InputStreamReader(System.in));
        pw = new PrintWriter(new BufferedWriter(new OutputStreamWriter(System.out)));
//        int s = readInt();
//        int [] ar = new int [s+1];
//        int [] ar1 = new int [s+1];
//        for (int i = 0; i <= s; i++){
//            ar[i] = readInt();
//        }
//        for (int i = 0; i <= s; i++){
//            ar1[i] = readInt();
//        }
//        Complex [] arr = pad(ar);
//        Complex [] arr1 = pad(ar1);
//        inPlace(arr, arr.length, false);
//        inPlace(arr1, arr1.length, false);
//        Complex [] mult = new Complex [arr.length];
//        for (int i = 0; i < arr.length; i++){
//            mult[i] = arr[i].multiply(arr1[i]);
//        }
//        for (int i = 0; i < mult.length; i++){
//            System.out.println(mult[i]);
//        }
//        inPlace(mult, mult.length, true);
//        for (int i = 0; i < s*2+1; i++){
//            long result = Math.round(mult[i].a/arr.length);
//            System.out.print(result + " ");
//
//        }
//        for (int s = 2; s <= 1000000*2; s <<=1) {
            int s = readInt();
            long[] a1 = new long[s + 1];
            long[] a2 = new long[s + 1];
            for (int i = 0; i <= s; i++) {
                a1[i] = readInt();
            }
            for (int i = 0; i <= s; i++) {
                a2[i] = readInt();
            }
            long[] ar = pad(a1);
            long[] ar1 = pad(a2);
            long[] ar2 = pad(a1);
            long[] ar3 = pad(a2);


//        int logSize = (int)Math.round(Math.log(ar.length)/Math.log(2));
//
//        long arr[] = findN(ar.length);
//        long arr2[] = findN2(ar.length, arr[1], arr[0]);
//
//        List<Long> primeFact1 = uniquePrimeFactors(arr[1] - 1);
//        List<Long> primeFact2 = uniquePrimeFactors(arr2[1] - 1);
//
//        BigInteger generator1 = findGenerator(primeFact1, arr[1]);
//        BigInteger generator2 = findGenerator(primeFact2, arr2[1]);
//
//        BigInteger precompute1 [] = new BigInteger[logSize];
//        BigInteger precompute2 [] = new BigInteger[logSize];
//        BigInteger precompute1I [] = new BigInteger[logSize];
//        BigInteger precompute2I [] = new BigInteger[logSize];
//
//        BigInteger kn1 = BigInteger.valueOf(arr[0]).multiply(BigInteger.valueOf(ar.length));
//        BigInteger kn2 = BigInteger.valueOf(arr2[0]).multiply(BigInteger.valueOf(ar.length));
//
////            System.out.println("Init time: " + (endTime-startTime)/1000000000.0);
//        int counter = 0;
//        for (int i = 2; i <= ar.length; i <<=1){
//            precompute1[counter] = findPrimitive(generator1,arr[1],kn1, i, false);
////                System.out.println(precompute1[counter]);
//            precompute2[counter] = findPrimitive(generator2, arr2[1], kn2, i, false);
//            precompute1I[counter] = findPrimitive(generator1,arr[1],kn1, i, true);
//            precompute2I[counter] = findPrimitive(generator2, arr2[1], kn2, i, true);
//            counter++;
//        }
//

//
//        BigInteger modInv = BigInteger.valueOf(arr2[1]).modInverse(BigInteger.valueOf(arr[1]));
////            BigInteger bigMod =

//        endTime = System.nanoTime();
////            System.out.println();
////            System.out.println((endTime-startTime)/1000000000.0 + "s");
////            System.out.println("---------------");
            long startTime = System.nanoTime();
            long root1[] = {2088763392, 78423909, 1596245477, 1571822688, 773715244, 641558202, 1603120306, 1876782378, 995984211, 583867446, 1522510314, 511837099, 1797322843, 859290800, 683478494, 845064300, 1956703310, 1993405836, 16814040, 49151182, 475662317, 12052012, 865472878};
            long root1I[] = {2088763392, 2010339484, 1503397660, 547032107, 1343641900, 447157216, 1769054244, 120157682, 667669475, 1423097324, 1489135927, 1724249164, 1108900612, 1104585648, 1985397649, 147791241, 129286670, 1217844064, 197092010, 792459346, 501317389, 1392927016, 1891284705};
            long root2[] = {2113929216, 911673634, 1550788396, 483614631, 324124180, 1852886862, 616074024, 152007057, 993942852, 505354142, 1402588562, 549981670, 2012629272, 1668209721, 1597367311, 477541733, 78760698, 1460160397, 1032064332, 263042423, 1976245982, 1331364166, 262581406};
            long root2I[] = {2113929216, 1202255583, 1824879522, 341521208, 1656523455, 2057797169, 791242072, 1995576413, 1061942124, 1420901561, 2081705832, 356482523, 940469005, 764304299, 1307976383, 1394223361, 746533687, 1832857487, 764557217, 1500383091, 400557881, 176344849, 882047861};
            long prime[] = {249, 2088763393};
            long prime2[] = {252, 2113929217};

            NTT(ar, ar.length, prime.clone(), root1.clone(),false);
            NTT(ar1, ar.length, prime.clone(), root1.clone(), false);
            long[] mult = new long[ar.length];
        for (int i = 0; i < ar.length; i++) {
//                mult[i] = BigInteger.valueOf(ar[i]).multiply(BigInteger.valueOf(ar1[i])).mod(BigInteger.valueOf(prime[1])).longValue();
            mult[i] = (ar[i]*ar1[i])%prime[1];
        }
        NTT(mult, ar.length, prime.clone(), root1I.clone(), true);
//        for (int i = 0; i < mult.length; i++) {
//            mult[i] = mult[i] / ar.length;
//        }
//        for (int i = 0; i < mult.length; i++){
//            System.out.println(mult[i] + " ");
//        }
            long[] mult1 = new long[ar.length];
            NTT(ar2, ar.length, prime2.clone(), root2.clone(), false);
            NTT(ar3, ar.length, prime2.clone(), root2.clone(), false);


            for (int i = 0; i < mult1.length; i++) {
//                mult1[i] = BigInteger.valueOf(ar2[i]).multiply(BigInteger.valueOf(ar3[i])).mod(BigInteger.valueOf(prime2[1])).longValue();
                mult1[i] = (ar2[i]*ar3[i])%prime2[1];
            }
            NTT(mult1, ar.length, prime2.clone(), root2I.clone(), true);
//        for (int i = 0; i < mult.length; i++) {
//            mult1[i] = mult1[i] / ar.length;
////        }
//        for (int i = 0; i < mult1.length; i++){
//            System.out.println(mult1[i]);
//        }
            BigInteger bigMod = BigInteger.valueOf(prime[1]).multiply(BigInteger.valueOf(prime2[1]));
            for (int i = 0; i < 2 * s + 1; i++) {
//            BigInteger diff = BigInteger.valueOf(mult[i]-mult1[i]);
//            if (diff.longValue() < 0){
//                System.out.println(diff);
//            }
//            BigInteger out = modInv.multiply(BigInteger.valueOf(prime[1])).multiply(diff).add(BigInteger.valueOf(mult[i]));
//                (crt(mult[i], mult1[i], prime[1], prime2[1], bigMod)).divide(BigInteger.valueOf(ar.length));
                System.out.print((crt(mult[i], mult1[i], prime[1], prime2[1], bigMod)) + " ");
//            if (out.longValue() < 0) {
//                System.out.print(out.longValue() + " ");
//            }

            }
            long endTime = System.nanoTime();
//            System.out.println(s);
//            System.out.println((endTime-startTime)/1000000000.0);
//        }

    }
}