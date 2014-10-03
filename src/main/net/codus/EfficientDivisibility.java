package net.codus;
import java.math.BigInteger;
import java.util.Arrays;


public class EfficientDivisibility implements Divisibility{
	
	private static final int BASE = 256*256;
	private static final OpenBigInteger BASE_BI = new OpenBigInteger("" + BASE);

	private OpenBigInteger[] multiplierByLastDigit = new OpenBigInteger[BASE];
	private int[] intMultiplierByLastDigit = new int[BASE];
	
	
	public EfficientDivisibility() {
		for (int i=1; i<BASE; i+=2){
			multiplierByLastDigit[i] = new OpenBigInteger(""+findMultiplyerByLastDigit(i, BASE));
		}
	}
	
//	private int getLastDigit(OpenBigInteger number){
//		byte[] numberBytes = number.toByteArray();
//		int lastDigit = (numberBytes[numberBytes.length-1] & 0x000000FF) + (numberBytes.length > 1 ? 256*(numberBytes[numberBytes.length-2] & 0x000000FF) : 0);
//		return lastDigit;
//	}
	
//	private OpenBigInteger divideByBase(OpenBigInteger number){
//		byte[] numberBytes = number.toByteArray();
//		if (numberBytes.length <=2 )
//			return OpenBigInteger.ZERO;
//		byte[] dividend = Arrays.copyOfRange(numberBytes, 0, numberBytes.length-2);
//		return new OpenBigInteger(dividend);
//	}
	
	@Override
	public boolean isDivisibile(BigInteger number, BigInteger prime) {
		return isDivisibile(new OpenBigInteger(number.toByteArray()), new OpenBigInteger(prime.toByteArray()));
	}
	
	// working now
	public boolean isDivisibile(OpenBigInteger number , OpenBigInteger prime){
		OpenBigInteger multiplyer = findMultiplyer(prime); 
		OpenBigInteger next = number;
		next.initializeMultiplyer(multiplyer);
		while (true){
			int lastDigit = next.getLastDigit();
//			OpenBigInteger numberLastDigit = OpenBigInteger.valueOf(lastDigit);
			next.divideByBase();
			next.subtractMultipication(multiplyer, lastDigit);
//			next = next.subtract(numberLastDigit.multiply(multiplyer));
//			System.out.println("Next :" + next.toString() +"  signum :" + next.signum());
			
			if (next.signum() <= 0 )
				break;
		}
		boolean result = next.toBigInteger().mod(prime.toBigInteger()).signum() == 0;
//		System.out.println("Result :" + result);
		return result;
	}
	
	private OpenBigInteger findMultiplyer(OpenBigInteger prime) {
		int lastDigit = prime.getLastDigit();
		OpenBigInteger multiplyer = multiplierByLastDigit[lastDigit];
		OpenBigInteger product = prime.multiply(multiplyer);
		product.divideByBase();
		return product;		
	}

//	Commented out since it requires the mod function
//	private static int findMultiplyerByTrial(OpenBigInteger prime) {
//		OpenBigInteger tmp = prime;
//		int mod = tmp.mod(BASE_BI).intValue();
////		System.out.println("prime :" + prime);
//		while (mod != 1){
////			System.out.println("mod :" + mod +" tmp :" + tmp);
//			tmp = tmp.add(prime).add(prime);
//			mod = tmp.mod(BASE_BI).intValue();
//		}
//		return tmp.divide(BASE_BI).intValue();
//	}
	
	private int findMultiplyerByLastDigit(int lastDigit, int base){
		long multiplyer = 1;
		long mod = lastDigit % base;
		long tmp = lastDigit;
		while (mod != 1){
			tmp += 2*lastDigit;
			multiplyer += 2;
			long lastMod = mod;
			mod = tmp % base;
			if (lastMod == mod)
				return -1;
		}
		if (multiplyer > Integer.MAX_VALUE)
			throw new RuntimeException("Multiplyer is not integer");
		return (int)multiplyer;
	}


	
//	public static void main(String[] args) {
//		int base = 256;
//		for (int i=1; i<base; i+=2){
//			System.out.println("i:" + i +":" + findMultiplyerByLastDigit(i, base));
//		}
//	}
	

//	public static void main(String[] args) {
//		isDivisibile(new OpenBigInteger("1642979"), new OpenBigInteger("47"));
//		isDivisibile(new OpenBigInteger("1642979"), new OpenBigInteger("47"));
//		isDivisibile(new OpenBigInteger("36249"), new OpenBigInteger("43"));
//		isDivisibile(new OpenBigInteger("925"), new OpenBigInteger("37"));
//		isDivisibile(new OpenBigInteger("837"), new OpenBigInteger("31"));
//		isDivisibile(new OpenBigInteger("14913708"), new OpenBigInteger("11"));
//		isDivisibile(new OpenBigInteger("867"), new OpenBigInteger("19"));
//	}

}
