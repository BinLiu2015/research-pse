package pseSimulation;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;

public class PSESimulator {
  private double[] allTaskInputs;
  private int taskNum = 400;
  private double stragglerPercent = 0.1; // how many tasks will become
                                         // straggler
  private boolean[] isSpeculative;
  private double B = 4; // copy bandwidth: MB/s
  private double selectivity = 0.05;
  private int slotPerNode = 4;
  private int nodeNum = 100;
  private int slotNum = slotPerNode * nodeNum;
  private double r = 0.02; // progressRate: will process I*r MB in one time unit
  private double[] slowPS; // when progress score reach this point,
                           // the chosen tasks run slowly;
  private int slowFactor = 10; // assume the straggling task will be 10 times
                               // slower

  public static double[] gaussian(int num, double mean, double dev) {
    Random R = new Random();
    double[] r = new double[num];
    for (int i = 0; i < num; i++) {
      r[i] = Math.sqrt(dev) * R.nextGaussian() + mean;
    }
    return r;
  }

  public void init() {
    taskNum = 400;
    allTaskInputs = gaussian(taskNum, 512, 64);
    stragglerPercent = 0.1;
    B = 4;
    selectivity = 0.05;
    slowFactor = 10;
    isSpeculative = new boolean[taskNum];
    slowPS = new double[taskNum];
    Random R = new Random();
    for (int i = 0; i < taskNum * stragglerPercent; i++) {
      int taskId = R.nextInt(taskNum);
      isSpeculative[taskId] = true;
      slowPS[taskId] = R.nextDouble();
      // System.out.println(taskId + ":" + slowPS[taskId]);
    }
    
    System.out.println("TaskNum: " + taskNum);
    System.out.println("All Task Inputs: " + Arrays.toString(allTaskInputs));
    System.out.println("Will be stragglers: " + Arrays.toString(isSpeculative));
    System.out.println("When start slowly: " + Arrays.toString(slowPS));
    System.out.println("Copy bandwidth: " + B);
    System.out.println("Slow Factor: " + slowFactor);
    System.out.println("Selectivity: " + selectivity);
  }

  // num of available slots = half of num of tasks
  public void smallCluster() {
    this.slotNum = taskNum % 2 == 0 ? taskNum / 2 : (taskNum + 1) / 2;

    int[] slots = new int[slotNum];
    int[] noSE = new int[taskNum];
    int[] late = new int[taskNum];
    int[] pse = new int[taskNum];
    int tid = 0;

    // NoSE:
    tid = 0;
    while (tid < slotNum) {
      int runTime = 0;
      if (isSpeculative[tid] == false) {
        runTime = ceilDouble(1 / r);
      } else {
        runTime = ceilDouble(slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
      }
      slots[tid] = runTime;
      noSE[tid] = runTime;
      tid++;
    }
    while (tid < taskNum) {
      // find a idle slot, launch a new task
      int runTime = 0;
      if (isSpeculative[tid] == false) {
        runTime = ceilDouble(1 / r);
      } else {
        runTime = ceilDouble(slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
      }
      noSE[tid] = runTime;
      int sid = findMinFromArray(slots);
      slots[sid] += runTime;
      tid++;
    }
    Arrays.sort(noSE);
    Arrays.sort(slots);
    System.out.println("NoSE:");
    System.out.println(Arrays.toString(noSE));
    System.out.println(Arrays.toString(slots));

    // LATE
    int lateSchedule = 0, lateSucc = 0;
    Map<Integer, Integer> slot2Task = new HashMap<Integer, Integer>();
    for (int i = 0; i < slotNum; i++) {
      slots[i] = 0;
    }
    tid = 0;
    while (tid < slotNum) {
      int runTime = 0;
      if (isSpeculative[tid] == false) {
        runTime = ceilDouble(1 / r);
      } else {
        runTime = ceilDouble(slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
      }
      slot2Task.put(tid, tid);
      slots[tid] = runTime;
      late[tid] = runTime;
      tid++;
    }
    int curTime = 0;
    while (tid < taskNum) {
      int runTime = 0;
      if (isSpeculative[tid] == false) {
        runTime = ceilDouble(1 / r);
      } else {
        runTime = ceilDouble(slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
      }
      late[tid] = runTime;
      int sid = findMinFromArray(slots);
      curTime = slots[sid];
      slots[sid] += runTime;
      slot2Task.put(sid, tid);
      tid++;
    }
    // slotId-->endTime
    Map<Integer, Integer> idleSlots = new HashMap<Integer, Integer>(); 
    // slotId-->endTime
    Map<Integer, Integer> busySlots = new HashMap<Integer, Integer>(); 
    // slotId -->endTime, need to be scheduled
    Map<Integer, Integer> slowSlots = new HashMap<Integer, Integer>(); 
    // slotId -->endTime, LATE estimated runtime
    Map<Integer, Integer> candidates = new HashMap<Integer, Integer>();
    Map<Integer, Integer> candidatesHasRun = new HashMap<Integer, Integer>();
    for (int sid = 0; sid < slotNum; sid++) {
      if (slots[sid] <= curTime) {
        idleSlots.put(sid, slots[sid]);
      } else {
        busySlots.put(sid, slots[sid]);
        if (isSpeculative[slot2Task.get(sid)]) {
          slowSlots.put(sid, slots[sid]);
          // LATE
          tid = slot2Task.get(sid);
          int shouldRun = ceilDouble(slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
          int phase1 = ceilDouble (slowPS[tid] / r);
          int realRun = curTime - (slots[sid] - shouldRun);
          candidates.put(sid, slots[sid] - curTime); // how long this slot/task
                                                     // will finish
          candidatesHasRun.put(sid, realRun);
        }
      }
    }
    while (!candidates.isEmpty()) {
      // should schedule
      int scheduledSlotId = findMinIdx(candidates);
      int idleSlotId = findMinIdx(idleSlots);
      assert idleSlots.get(idleSlotId) <= curTime;  //assuming there is enough idle slots

      int oriEndTime = slots[scheduledSlotId];
      int speEndTime = curTime + ceilDouble(1 / r);
      lateSchedule++;
      int endTime = oriEndTime;
      if (speEndTime < oriEndTime) {
        endTime = speEndTime;
        lateSucc++;
        //the speculative succeeds
        late[slot2Task.get(scheduledSlotId)] = candidatesHasRun.get(scheduledSlotId) + ceilDouble (1 / r);  
      }
      idleSlots.put(idleSlotId, endTime);
      slots[scheduledSlotId] = endTime;
      slots[idleSlotId] = endTime;   //this slot may be idle for a while between two tasks;
      candidates.remove(scheduledSlotId);
    }

     Arrays.sort(late);
     Arrays.sort(slots);
     System.out.println("LATE: " + lateSucc + "/" + lateSchedule);
     System.out.println(Arrays.toString(late));
     System.out.println(Arrays.toString(slots));
     
     // PSE
     int pseSchedule = 0, pseSucc = 0;
     slot2Task.clear();
     for (int i = 0; i < slotNum; i++) {
       slots[i] = 0;
     }
     tid = 0;
     while (tid < slotNum) {
       int runTime = 0;
       if (isSpeculative[tid] == false) {
         runTime = ceilDouble(1 / r);
       } else {
         runTime = ceilDouble(slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
       }
       slot2Task.put(tid, tid);
       slots[tid] = runTime;
       pse[tid] = runTime;
       tid++;
     }
     curTime = 0;
     while (tid < taskNum) {
       int runTime = 0;
       if (isSpeculative[tid] == false) {
         runTime = ceilDouble(1 / r);
       } else {
         runTime = ceilDouble(slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
       }
       pse[tid] = runTime;
       int sid = findMinFromArray(slots);
       curTime = slots[sid];
       slots[sid] += runTime;
       slot2Task.put(sid, tid);
       tid++;
     }
     // slotId-->endTime
     idleSlots.clear();
     busySlots.clear();
     slowSlots.clear();
     candidates.clear();
     //How long these candidates will finish since curTime
     Map<Integer, Integer> candidatesRT = new HashMap<Integer, Integer>();
     for (int sid = 0; sid < slotNum; sid++) {
       if (slots[sid] <= curTime) {
         idleSlots.put(sid, slots[sid]);
       } else {
         busySlots.put(sid, slots[sid]);
         if (isSpeculative[slot2Task.get(sid)]) {
           slowSlots.put(sid, slots[sid]);
           // PSE
           int oriRT = slots[sid] - curTime;
           tid = slot2Task.get(sid);
           int shouldRun = ceilDouble(slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
           int phase1 = ceilDouble (slowPS[tid] / r);
           int realRun = curTime - (slots[sid] - shouldRun);
           double processedInputData = 0;
           if(realRun <= phase1) {
             processedInputData = realRun * r;
           } else {
             processedInputData = allTaskInputs[tid]*slowPS[tid] + (realRun-phase1)*r/slowFactor;
           }
           int pseRT =
               ceilDouble (processedInputData * selectivity / B + ((allTaskInputs[tid] - processedInputData) / allTaskInputs[tid])
                  / r);
          if (oriRT - pseRT > 0) {
            candidates.put(sid, oriRT - pseRT); // speculation gain
            candidatesRT.put(sid, pseRT);
            pse[tid] = realRun + pseRT;  //Accutual run time for this task;
          }
         }
       }
     }
     while (!candidates.isEmpty()) {
       // should schedule
       int scheduledSlotId = findMaxIdx(candidates);  //maximum gain
       int idleSlotId = findMinIdx(idleSlots);
       assert idleSlots.get(idleSlotId) <= curTime;  //assuming there is enough idle slots

       pseSchedule++;
       pseSucc++;
       int endTime = curTime + candidatesRT.get(scheduledSlotId);
       //pse[slot2Task.get(scheduledSlotId)] has been assigned when estimating gain
       idleSlots.put(idleSlotId, endTime);
       slots[scheduledSlotId] = endTime;
       slots[idleSlotId] = endTime;   //this slot may be idle for a while between two tasks;
       candidates.remove(scheduledSlotId);
     }
     Arrays.sort(pse);
     Arrays.sort(slots);
     System.out.println("PSE: " + pseSucc + "/" + pseSchedule);
     System.out.println(Arrays.toString(pse));
     System.out.println(Arrays.toString(slots));
  }

  // num of available nodes = num of tasks
  public void mediumCluster() {
    int[] noSE = new int[taskNum];
    int[] se = new int[taskNum];
    int[] late = new int[taskNum];
    int[] pse = new int[taskNum];

    // NoSE
    for (int tid = 0; tid < taskNum; tid++) {
      if (isSpeculative[tid] == false) {
        noSE[tid] = ceilDouble (1 / r);
      } else {
        noSE[tid] = ceilDouble(slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
      }
    }

    // TODO: SE

    // LATE
    int lateSchedule = 0, lateSucc = 0;
    Map<Integer, Integer> normalTasks = new HashMap<Integer, Integer>();
    Map<Integer, Integer> slowTasks = new HashMap<Integer, Integer>();
    Map<Integer, Double> slowTaskProgressScore = new HashMap<Integer, Double>();
    for (int tid = 0; tid < taskNum; tid++) {
      if (isSpeculative[tid] == false) {
        late[tid] = ceilDouble (1 / r);
        normalTasks.put(tid, late[tid]);
      } else {
        late[tid] = ceilDouble (slowPS[tid] / r);
        slowTasks.put(tid, late[tid]);
        slowTaskProgressScore.put(tid, slowPS[tid]);
      }
    }
    Map<Integer, Double> candidates = new HashMap<Integer, Double>();
    while (!slowTasks.isEmpty()) {
      // System.out.println(slowTasks);
      // System.out.println(Arrays.toString(slowTasks.values().toArray()));
      candidates.clear();
      int id1 = findMinIdx(normalTasks);
      int id2 = findMinIdx(slowTasks);
      int t1 = normalTasks.get(id1); // the earliest idle slot
      int t2 = slowTasks.get(id2);
      if (t2 < t1) {
        // no idle slot: the slow tasks whose time is earlier than t1 have to
        // wait another (t1-t2)
        for (Entry<Integer, Integer> e : slowTasks.entrySet()) {
          if (e.getValue() < t1) {// compute candidates
            double addedPS = ((t1 - e.getValue()) * r / slowFactor) /*
                                                                     * /
                                                                     * allTaskInputs
                                                                     * [
                                                                     * e.getKey(
                                                                     * )]
                                                                     */;
            candidates.put(e.getKey(), slowTaskProgressScore.get(e.getKey())
                + addedPS);
          }
        }
        // So far, can schedule at least one
        int scheduledTask = findMinIdx_2(candidates);
        late[scheduledTask] = t1;
        int oriNeed = ceilDouble ((1 - candidates.get(scheduledTask)) / (r / slowFactor));
        int speNeed = ceilDouble (1 / r);
        int needTime = oriNeed;
        lateSchedule++;
        if (needTime > speNeed) {
          needTime = speNeed;
          lateSucc++;
        }
        slowTasks.remove(scheduledTask); // schedule it
        slowTaskProgressScore.remove(scheduledTask);
        late[scheduledTask] += needTime; // Note that only the scheduled one
                                         // should be updated!
        normalTasks.put(id1, t1 + needTime); // when it will be idle
      } else { // t1 <= t2
        // can schedule one right now
        for (Entry<Integer, Integer> e : slowTasks.entrySet()) {
          if (e.getValue() == t2) {
            candidates.put(e.getKey(), slowTaskProgressScore.get(e.getKey()));
          }
        }
        int scheduledTask = findMinIdx_2(candidates);
        int oriNeed = ceilDouble ((1 - candidates.get(scheduledTask)) / (r / slowFactor));
        int speNeed = ceilDouble (1 / r) ;
        int needTime = oriNeed;
        lateSchedule++;
        if (needTime > speNeed) {
          needTime = speNeed;
          lateSucc++;
        }
        slowTasks.remove(scheduledTask); // schedule it
        slowTaskProgressScore.remove(scheduledTask);
        late[scheduledTask] += needTime;
        normalTasks.put(id1, t2 + needTime); // when it will be idle
      }
    }

    // PSE
    int pseSchedule = 0, pseSucc = 0;
    normalTasks.clear();
    slowTasks.clear();
    slowTaskProgressScore.clear();
    for (int tid = 0; tid < taskNum; tid++) {
      if (isSpeculative[tid] == false) {
        pse[tid] = ceilDouble(1 / r);
        normalTasks.put(tid, pse[tid]);
      } else {
        pse[tid] = ceilDouble(slowPS[tid] / r);
        slowTasks.put(tid, pse[tid]);
        slowTaskProgressScore.put(tid, slowPS[tid]);
      }
    }
    Map<Integer, Double> candidatesPS = new HashMap<Integer, Double>();
    Map<Integer, Double> candidatesGain = new HashMap<Integer, Double>();
    while (!slowTasks.isEmpty()) {
      // System.out.println(slowTasks);
      // System.out.println(Arrays.toString(slowTasks.values().toArray()));
      candidatesPS.clear();
      candidatesGain.clear();
      int id1 = findMinIdx(normalTasks);
      int id2 = findMinIdx(slowTasks);
      int t1 = normalTasks.get(id1); // the earliest idle slot
      int t2 = slowTasks.get(id2);
      if (t2 < t1) {
        // no idle slot: the slow tasks whose time is earlier than t1 have to
        // wait another (t1-t2)
        for (Entry<Integer, Integer> e : slowTasks.entrySet()) {
          if (e.getValue() < t1) {// compute candidates
            double addedPS = ((t1 - e.getValue()) * r / slowFactor) /*
                                                                     * /
                                                                     * allTaskInputs
                                                                     * [
                                                                     * e.getKey(
                                                                     * )]
                                                                     */;
            candidatesPS.put(e.getKey(), slowTaskProgressScore.get(e.getKey())
                + addedPS);
            // speculation gain = t_ori - t_pse
            double t_ori =
                (1 - candidatesPS.get(e.getKey())) / (r / slowFactor);
            double t_pse =
                allTaskInputs[e.getKey()] * candidatesPS.get(e.getKey())
                    * selectivity / B + (1 - candidatesPS.get(e.getKey())) / r;
            double gain = t_ori - t_pse;
            candidatesGain.put(e.getKey(), gain);
          }
        }
        // So far, schedule one if there is a candidate
        int scheduledTask = findMaxIdx_2(candidatesGain);
        if (scheduledTask < 0) {
          // no candidate
          for (Entry<Integer, Integer> e : slowTasks.entrySet()) {
            if (e.getValue() < t1) { // update progress score
              double addedPS = ((t1 - e.getValue()) * r / slowFactor) /*
                                                                       * /
                                                                       * allTaskInputs
                                                                       * [
                                                                       * e.getKey
                                                                       * ()]
                                                                       */;
              pse[e.getKey()] += (t1 - e.getValue());
              slowTasks.put(e.getKey(), pse[e.getKey()]);
              slowTaskProgressScore.put(e.getKey(),
                  slowTaskProgressScore.get(e.getKey()) + addedPS);
            }
          }
          continue;
        }
        pse[scheduledTask] = t1;
        int oriNeed =ceilDouble ((1 - candidatesPS.get(scheduledTask)) / (r / slowFactor));
        int speNeed =ceilDouble (allTaskInputs[scheduledTask]
                * candidatesPS.get(scheduledTask) * selectivity / B + (1 - candidatesPS
                .get(scheduledTask)) / r);
        int needTime = oriNeed;
        pseSchedule++;
        if (needTime > speNeed) {
          needTime = speNeed;
          pseSucc++;
        }
        slowTasks.remove(scheduledTask); // schedule it
        slowTaskProgressScore.remove(scheduledTask);
        pse[scheduledTask] += needTime; // Note that only the scheduled one
                                        // should be updated!
        normalTasks.put(id1, t1 + needTime); // when it will be idle
      } else { // t1 <= t2
        // can schedule one right now
        for (Entry<Integer, Integer> e : slowTasks.entrySet()) {
          if (e.getValue() == t2) {
            candidatesPS.put(e.getKey(), slowTaskProgressScore.get(e.getKey()));
            // speculation gain = t_ori - t_pse
            double t_ori =
                (1 - candidatesPS.get(e.getKey())) / (r / slowFactor);
            double t_pse =
                allTaskInputs[e.getKey()] * candidatesPS.get(e.getKey())
                    * selectivity / B + (1 - candidatesPS.get(e.getKey())) / r;
            double gain = t_ori - t_pse;
            candidatesGain.put(e.getKey(), gain);
          }
        }
        int scheduledTask = findMaxIdx_2(candidatesGain);
        if (scheduledTask < 0) {
          // no candidate: remove id2, it will use original task
          pse[id2] += ceilDouble((1 - slowTaskProgressScore.get(id2)) / (r / slowFactor));
          slowTasks.remove(id2);
          slowTaskProgressScore.remove(id2);
          continue;
        }
        int oriNeed = ceilDouble((1 - candidates.get(scheduledTask)) / (r / slowFactor));
        int speNeed = ceilDouble(allTaskInputs[scheduledTask]
                * candidatesPS.get(scheduledTask) * selectivity / B + (1 - candidatesPS
                .get(scheduledTask)) / r);
        int needTime = oriNeed;
        pseSchedule++;
        if (needTime > speNeed) {
          needTime = speNeed;
          pseSucc++;
        }
        slowTasks.remove(scheduledTask); // schedule it
        slowTaskProgressScore.remove(scheduledTask);
        pse[scheduledTask] += needTime; // Note that only the scheduled one
                                        // should be updated!
        normalTasks.put(id1, t1 + needTime); // when it will be idle
      }
    }

    Arrays.sort(noSE);
    Arrays.sort(late);
    Arrays.sort(pse);

    System.out.println("NoSE:");
    System.out.println(Arrays.toString(noSE));

    System.out.println("LATE: " + lateSucc + "/" + lateSchedule);
    System.out.println(Arrays.toString(late));

    System.out.println("PSE: " + pseSucc + "/" + pseSchedule);
    System.out.println(Arrays.toString(pse));
  }

  // num of available nodes >= num of tasks + num of slow tasks
  public void largeCluster() {
    int[] noSE = new int[taskNum];
    int[] se = new int[taskNum];
    int[] late = new int[taskNum];
    int[] pse = new int[taskNum];

    int lateSchedule = 0, lateSucc = 0;
    int pseSchedule = 0, pseSucc = 0;

    for (int tid = 0; tid < taskNum; tid++) {
      if (isSpeculative[tid] == false) {
        noSE[tid] = ceilDouble(1 / r);
        se[tid] = ceilDouble(1 / r);
        late[tid] = ceilDouble(1 / r);
        pse[tid] = ceilDouble(1 / r);
      } else {
        noSE[tid] = ceilDouble (slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
        se[tid] = ceilDouble (slowPS[tid] / r + 1 / r);

        double ori = (1 - slowPS[tid]) / (r / slowFactor);
        double full = 1 / r;
        lateSchedule++;
        if (ori < full) {
          late[tid] =ceilDouble(slowPS[tid] / r + (1 - slowPS[tid]) / (r / slowFactor));
        } else {
          late[tid] = ceilDouble(slowPS[tid] / r + 1 / r);
          lateSucc++;
        }

        double partial =
            allTaskInputs[tid] * slowPS[tid] * selectivity / B
                + (1 - slowPS[tid]) / r;
        double needTime = full;
        if (needTime > partial) {
          needTime = partial;
          pseSchedule++;
          pseSucc++;
        }
        pse[tid] = ceilDouble(slowPS[tid] / r + needTime);
      }
    }
    Arrays.sort(noSE);
    Arrays.sort(late);
    Arrays.sort(pse);

    System.out.println("NoSE:");
    System.out.println(Arrays.toString(noSE));

    System.out.println("LATE: " + lateSucc + "/" + lateSchedule);
    System.out.println(Arrays.toString(late));

    System.out.println("LATE: " + pseSucc + "/" + pseSchedule);
    System.out.println(Arrays.toString(pse));
  }

  private int findMinFromArray(int[] arr) {
    int idx = -1, minV = Integer.MAX_VALUE;
    for (int i = 0; i < arr.length; i++) {
      if (arr[i] < minV) {
        idx = i;
        minV = arr[i];
      }
    }
    return idx;
  }

  private int findMinIdx(Map<Integer, Integer> id2Time) {
    Iterator<Integer> ite = id2Time.keySet().iterator();
    int idx = -1, minV = Integer.MAX_VALUE;
    while (ite.hasNext()) {
      int id = ite.next();
      if (id2Time.get(id) < minV) {
        idx = id;
        minV = id2Time.get(idx);
      }
    }
    return idx;
  }

  private int findMaxIdx(Map<Integer, Integer> gain) {
    Iterator<Integer> ite = gain.keySet().iterator();
    int idx = -1, maxV = Integer.MIN_VALUE;
    while (ite.hasNext()) {
      int id = ite.next();
      if (gain.get(id) > maxV) {
        idx = id;
        maxV = gain.get(idx);
      }
    }
    return idx;
  }
  
  private int findMinIdx_2(Map<Integer, Double> id2Time) {
    Iterator<Integer> ite = id2Time.keySet().iterator();
    int idx = -1;
    double minV = Double.MAX_VALUE;
    while (ite.hasNext()) {
      int id = ite.next();
      if (id2Time.get(id) < minV) {
        idx = id;
        minV = id2Time.get(idx);
      }
    }
    return idx;
  }

  private int findMaxIdx_2(Map<Integer, Double> id2Gain) {
    Iterator<Integer> ite = id2Gain.keySet().iterator();
    int idx = -1;
    double maxV = Double.MIN_VALUE;
    while (ite.hasNext()) {
      int id = ite.next();
      if (id2Gain.get(id) > maxV) {
        idx = id;
        maxV = id2Gain.get(idx);
      }
    }
    return idx;
  }

  public static int ceilDouble(double v) {
    int r = (int) v;
    double diff = v - r;
    if(diff > 0.00000000001){
      r++;
    }
    return r;
  }

  public static void main(String[] args) {
    // double[] r = gaussian(20, 10, 1);
    // System.out.println(Arrays.toString(r));
//    System.out.println(ceilDouble(1.0));
    
    PSESimulator obj = new PSESimulator();
    obj.init();

    System.out
        .println("\nSmall Cluster: ==============================================");
    obj.smallCluster();

    System.out
        .println("\nMedium Cluster: ==============================================");
    obj.mediumCluster();

    System.out
        .println("\nLarge Cluster: ==============================================");
    obj.largeCluster();
  }

}
