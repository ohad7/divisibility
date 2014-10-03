package net.codus;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;


public class BaseXDivisiblity implements Divisibility{

	private static final int BASE = 256*256;
	private static final BigInteger BASE_BI = new BigInteger("" + BASE);

	private BigInteger[] multiplierByLastDigit = new BigInteger[BASE];
	
	
	public BaseXDivisiblity() {
		for (int i=1; i<BASE; i+=2){
			multiplierByLastDigit[i] = new BigInteger(""+findMultiplyerByLastDigit(i, BASE));
		}
	}
	
	@Override
	public boolean isDivisibile(BigInteger number , BigInteger prime){
		BigInteger multiplyer = findMultiplyer(prime); 
		BigInteger next = number;
		while (true){
			BigInteger numberLastDigit = next.mod(BASE_BI);
			next = next.divide(BASE_BI).subtract(numberLastDigit.multiply(multiplyer));
//			System.out.println("Next :" + next.toString() +"  signum :" + next.signum());
			
			if (next.signum() <= 0 )
				break;
		}
		boolean result = next.mod(prime).intValue() == 0;
//		System.out.println("Result :" + result);
		return result;
	}
	
	private BigInteger findMultiplyer(BigInteger prime) {
		int lastDigit = prime.mod(BASE_BI).intValue();
		BigInteger multiplyer = multiplierByLastDigit[lastDigit];
		return prime.multiply(multiplyer).divide(BASE_BI);
	}

	private static int findMultiplyerByTrial(BigInteger prime) {
		BigInteger tmp = prime;
		int mod = tmp.mod(BASE_BI).intValue();
//		System.out.println("prime :" + prime);
		while (mod != 1){
//			System.out.println("mod :" + mod +" tmp :" + tmp);
			tmp = tmp.add(prime).add(prime);
			mod = tmp.mod(BASE_BI).intValue();
		}
		return tmp.divide(BASE_BI).intValue();
	}
	
	private int findMultiplyerByLastDigit(int lastDigit, int base){
		int multiplyer = 1;
		int mod = lastDigit % base;
		int tmp = lastDigit;
		while (mod != 1){
			tmp += 2*lastDigit;
			multiplyer += 2;
			int lastMod = mod;
			mod = tmp % base;
			if (lastMod == mod)
				return -1;
		}
		return multiplyer;
	}
	
//	public static void main(String[] args) {
//		int base = 256;
//		for (int i=1; i<base; i+=2){
//			System.out.println("i:" + i +":" + findMultiplyerByLastDigit(i, base));
//		}
//	}
	

//	public static void main(String[] args) {
//		isDivisibile(new BigInteger("1642979"), new BigInteger("47"));
//		isDivisibile(new BigInteger("1642979"), new BigInteger("47"));
//		isDivisibile(new BigInteger("36249"), new BigInteger("43"));
//		isDivisibile(new BigInteger("925"), new BigInteger("37"));
//		isDivisibile(new BigInteger("837"), new BigInteger("31"));
//		isDivisibile(new BigInteger("14913708"), new BigInteger("11"));
//		isDivisibile(new BigInteger("867"), new BigInteger("19"));
//	}
	
}
